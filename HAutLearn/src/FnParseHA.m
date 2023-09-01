function conditions = FnParseHA(xmlstruct, ode, label_guard)
% FnParseHA converts the xml representation of an HA created by HAutLearn
% to a condition based representation. State ids are changed to align with
% the the ode ids.

    % Create mapping between location ids from xml-representation and ode ids
    locToClusterID = computeMapLocToCluster(xmlstruct,ode);
    
    % Extract conditions from xml-representation and transform ids
    conditions = computeConditions(xmlstruct,label_guard,locToClusterID);
end

function locToClusterID = computeMapLocToCluster(xmlstruct,ode)
% computeMapLocToCluster generates a mapping from location ids
% (xml-representation) to ode ids
global num_var num_ud offsetCluster
    
    locations = xmlstruct.component(1).location;
    locToClusterID = containers.Map('KeyType','double','ValueType','double');

    % Find for each location corresponding entry in ode
    for i = 1:length(locations)
        % Convert string representation for currently considered location to a system matrix 
        flow = locations(i).flow;
        num_vars = num_var * (1+offsetCluster);
        AB = zeros(num_vars,num_vars+num_ud+1);
        lines = split(flow,'&');

        % Each line of string representation contributes one row to system matrix
        for j = 1:num_vars
            curr_line = lines(j);

            % Remove certain spaces to associate sign with coefficient 
            curr_line = strrep(curr_line,"+ ","+");
            curr_line = char(strrep(curr_line,"- ","-"));

            % Extract coefficients for all considered variables
            [AB,curr_line] = extractCoeffs(AB,curr_line,j,"x",num_vars,0);
            [AB,curr_line] = extractCoeffs(AB,curr_line,j,"u",num_ud,num_vars);

            % Extract constant in equation
            pos_end = length(curr_line);
            pos_start = strfind(curr_line(1:(pos_end-1))," ");
            pos_start = pos_start(end);
            coeff = sscanf(curr_line((pos_start):pos_end),"%f");
            % No constant present
            if(isempty(coeff))
                continue;
            end
            % Save constant in system matrix
            AB(j,num_vars+num_ud+1) = coeff;
        end

        % Match extracted system matrix to closest entry in ode
        dev = zeros(1,length(ode));
        for j = 1:length(ode)
            dev(1,j) = sum(sum(abs(cell2mat(ode(1,j))-AB))); %maybe other matrix norm
        end
        % There is exactly one match represented by minimal deviation
        [mindev,pos_min] = min(dev);
        % Store in map: key location id (i) and value ode id (j)
        locToClusterID(i) = pos_min;
    end
end

function [AB,curr_line] = extractCoeffs(AB,curr_line,j,c,num,off)
% extractCoeffs extracts all coeffs of form "character c followed k" (up
% to num) from curr_line and saves them in row j and col off+k of AB
    for k = 1:num
        pos_end = strfind(curr_line,sprintf(" * "+c+"%d",k))-1;
        % Considered variable not present, thus no coefficient needed
        if(isempty(pos_end))
            continue;
        end
        % Considered variable present, thus extract coefficient
        pos_start = strfind(curr_line(1:pos_end)," ");
        pos_start = pos_start(end);
        coeff = sscanf(curr_line((pos_start):pos_end),"%f");
        curr_line((pos_start):(pos_end+5)) = [];
        if(isempty(coeff))
            continue;
        end
        % Save coefficient in respective entry in sytem matrix
        AB(j,off+k) = coeff;
    end
end

function conditions = computeConditions(xmlstruct,label_guard,locToClusterID)
% computeConditions converts all transitions from an xml-representation into
% a condition based representation with ids aligning with the ode numbering
global num_var num_ud
    
    transitions = xmlstruct.component(1).transition;
    num = 1;
    conditions = [];
    
    % Consider all transitions of the xml-representation
    for i = 1:length(xmlstruct.component(1).transition)
        curr_trans = transitions(i); 
        condition = zeros(1,2+num_var+num_ud+1);

        % Extract Origin
        condition(1,1) = locToClusterID(curr_trans.sourceAttribute);
        % Extract Destination
        condition(1,2) = locToClusterID(curr_trans.targetAttribute);
        % Remove self transitions
        if (condition(1,1) == condition(1,2))
            continue;
        end

        curr_line = curr_trans.guard;
        % Remove certain spaces to associate sign with coefficient
        curr_line = strrep(curr_line,"+ ","+");
        curr_line = char(strrep(curr_line,"- ","-"));
        curr_line = [' ' curr_line];

        % Extract coefficients for all considered variables
        [condition,curr_line] = extractCoeffs(condition,curr_line,1,"x",num_var,2);
        [condition,curr_line] = extractCoeffs(condition,curr_line,1,"u",num_ud,2+num_var);

        % Extract remaining part of condition
        if(contains(curr_line,">"))
            condition(1,2+num_var+num_ud+1) = +1; %representing +1.0 > 0.0
        else
            condition(1,2+num_var+num_ud+1) = -1; %representing +1.0 < 0.0
        end

        % Match extracted condition vector to closest entry in label_guard
        dev = zeros(1,length(label_guard));
        for j = 1:length(label_guard)
            % The following comment is wrong, CONTINUE HERE AARON:
            % Constants is always +1.0 in normal form, thus not used in the comparison 
            comp_cond = condition(1,3:end);
            comp_guard = cell2mat(label_guard(1,j))';
            comp_guard(num_var+num_ud+1) = comp_guard(end);
            comp_guard = comp_guard(1:length(comp_cond));
            % Compute deviation between current condition vector and current label guard
            dev(1,j) = sum(sum(abs(comp_cond-comp_guard))); %maybe other matrix norm
        end
        % There is exactly one match represented by minimal deviation
        [mindev,pos_min] = min(dev);

        % Enhance condition with unrounded entry out of label_guard
        enhance_guard = cell2mat(label_guard(1,pos_min))';
        condition(1,3:(num_var+num_ud+2)) = enhance_guard(1,1:(num_var+num_ud));
        conditions = [conditions; condition];
        num = num + 1;
    end
end