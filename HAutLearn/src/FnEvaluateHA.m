function [correct,false,conditions] = FnEvaluateHA(trace, xmlstruct, ode, tol, label_guard)
global num_var offsetCluster
    %compute relation between clustering ids and location ids
    locations = xmlstruct.component(1).location;
    locToClusterID = containers.Map('KeyType','double','ValueType','double');
    
    for i = 1:length(locations)
        flow = locations(i).flow;
        num_vars = num_var * (1+offsetCluster);
        AB = zeros(num_vars,num_vars+1);
        lines = split(flow,'&');
        %create system matrix for currently considered location
        for j = 1:num_vars
            curr_line = lines(j);
            curr_line = strrep(curr_line,"+ ","+");
            curr_line = char(strrep(curr_line,"- ","-"));
            for k = 1:num_vars
                pos_end = strfind(curr_line,sprintf(" * x%d",k))-1;
                if(isempty(pos_end))
                    continue;
                end
                pos_start = strfind(curr_line(1:pos_end)," ");
                pos_start = pos_start(end);
                coeff = sscanf(curr_line((pos_start):pos_end),"%f");
                curr_line((pos_start):(pos_end+5)) = [];
                if(isempty(coeff))
                    continue;
                end
                AB(j,k) = coeff;
            end
            pos_end = length(curr_line);
            pos_start = strfind(curr_line(1:(pos_end-1))," ");
            pos_start = pos_start(end);
            coeff = sscanf(curr_line((pos_start):pos_end),"%f");
            if(isempty(coeff))
                continue;
            end
            AB(j,num_vars+1) = coeff;
        end
        %find corresponding matrix in ode (within tolerance)
        dev = zeros(1,length(ode));
        for j = 1:length(ode)
            dev(1,j) = sum(sum(abs(cell2mat(ode(1,j))-AB))); %maybe other matrix norm
        end
        %store in map, key location id (i) and value ode id (j)
        [mindev,pos_min] = min(dev);
        locToClusterID(i) = pos_min;
    end
    
    %extract transitions
    transitions = xmlstruct.component(1).transition;
    num = 1;
    conditions = [];
    
    for i = 1:length(xmlstruct.component(1).transition)
        curr_trans = transitions(i); 
        condition = zeros(1,num_var+1+2);
        condition(1,1) = locToClusterID(curr_trans.sourceAttribute); %origin
        condition(1,2) = locToClusterID(curr_trans.targetAttribute); %destination
        %remove self transitions
        if (condition(1,1) == condition(1,2))
            continue;
        end
        curr_line = curr_trans.guard;
        curr_line = strrep(curr_line,"+ ","+");
        curr_line = char(strrep(curr_line,"- ","-"));
        curr_line = [' ' curr_line];
        for k = 1:num_var
            pos_end = strfind(curr_line,sprintf(" * x%d",k))-1;
            if(isempty(pos_end))
                continue
            end
            pos_start = strfind(curr_line(1:pos_end)," ");
            pos_start = pos_start(end);
            coeff = sscanf(curr_line((pos_start):pos_end),"%f");
            curr_line((pos_start):(pos_end+5)) = [];
            condition(1,k+2) = coeff;
        end
        if(contains(curr_line,">"))
            condition(1,2+num_var+1) = +1; %representing +1.0 > 0.0
        else
            condition(1,2+num_var+1) = -1; %representing +1.0 < 0.0
        end
        %find corresponding vector in label_guard (within tolerance)
        dev = zeros(1,length(label_guard));
        for j = 1:length(label_guard)
            comp_cond = condition(1,3:end);
            comp_guard = cell2mat(label_guard(1,j))';
            %trim constants out of label guard
            comp_guard(num_var+1) = comp_guard(end);
            comp_guard = comp_guard(1:length(comp_cond));
            dev(1,j) = sum(sum(abs(comp_cond-comp_guard))); %maybe other matrix norm
        end
        %enhance condition with unrounded entry out of label_guard
        [mindev,pos_min] = min(dev);
        enhance_guard = cell2mat(label_guard(1,pos_min))';
        condition(1,3:(num_var+2)) = enhance_guard(1,1:num_var);
        conditions = [conditions; condition];
        num = num + 1;
    end
    
    %evaluate results
    correct = 0;
    false = 0;
    %HAutLearn use CP as last entry of segment, while we use it as first entry
    %of new segment, thus HAutLearn Notation used here (maybe change own
    %notation for better comparisons)
    for l = 1:length(trace)
        state_index = 1;
        for i = 1:(size(trace(l).x,1)-1)
            %last changepoints last datapoint thus no boundary check needed
            if(i>trace(l).chpoints(state_index+1)) 
                state_index = state_index + 1;
            end
            state = trace(l).labels_trace(state_index);
            predicted_state = state;
            for j = 1:size(conditions,1) %determinsm and no multiple trans assumed
                condition = conditions(j,:);
                if(condition(1,1) == state)
                    overall_sum = +1.0; %always part of condition generated by HAutLearn
                    for k = 1:num_var
                        overall_sum = overall_sum + trace(l).x(i,k) * condition(1,2+k);
                    end
                    %condition holds
                    if((overall_sum < tol && condition(1,2+num_var+1) == -1) || ...
                       (overall_sum > -tol && condition(1,2+num_var+1) == +1))
                        predicted_state = condition(1,2);
                        break;
                    end
                end
            end
            %compare to actual next state
            next_state = state;
            if(i+1>trace(l).chpoints(state_index+1))
                next_state = trace(l).labels_trace(state_index+1);
            end
            %predicition correct
            if(next_state == predicted_state)
                correct = correct + 1;
            else
                false = false + 1;
                %disp([num2str(overall_sum),int2str(condition(1,2+num_var+1))])
            end
        end
    end
end