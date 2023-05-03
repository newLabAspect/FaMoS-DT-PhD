function ode = FnEstODE(trace)
global Ts num_var num_ud offsetCluster
len_labels = length(trace(1).labels_num);
for label = 1:len_labels
    max_num_var = num_var * (1 + offsetCluster);
    for off = offsetCluster:-1:0
        x_seg = [];
        x_seg_plus = [];
        ud_seg = [];
        num_vars = num_var * (1+off);
        for j = 1:length(trace)
            labels_trace = trace(j).labels_trace;
            idx = find(labels_trace == label);
            x = trace(j).x(:,1:num_vars);
            ud = trace(j).ud;
             
            startj = trace(j).chpoints(idx);
            endj = trace(j).chpoints(idx+1)-1;
            
            for n = 1:length(startj)
                x_seg = [x_seg, x(startj(n):(endj(n)-1), :)'];
                if isempty(ud)
                    ud_seg = [ud_seg, ones(1, (endj(n)- startj(n)))];
                else
                    ud_seg = [ud_seg, [ud(startj(n):(endj(n)-1), :)'; ones(1, (endj(n)- startj(n)))]];
                end
                x_seg_plus = [x_seg_plus, x((startj(n)+1):endj(n), :)'];
            end
        end
        
        if size(x_seg,2)/length(trace)==2
            %b is not estimated
            if num_ud == 0
                A_Bu = mrdivide(x_seg_plus,[x_seg]);
                dA = A_Bu(:,1:num_vars);
                dB = zeros(size(A_Bu,1),1);
                % -1 indicate the ode is descrete and is used for variable reset
                ode(label) = {[padarray(dA,[max_num_var-num_vars max_num_var-num_vars],0,'post'), padarray(dB,[max_num_var-num_vars],0,'post'), -ones(num_vars,1)]};
            %b is estimated
            else
                A_Bu = mrdivide(x_seg_plus,[x_seg;ud_seg]);
                dA = A_Bu(:,1:num_vars);
                dB = A_Bu(:,num_vars+1:end);
                % -1 indicate the ode is descrete and is used for variable reset
                ode(label) = {[padarray(dA,[max_num_var-num_vars max_num_var-num_vars],0,'post'), padarray(dB,[max_num_var-num_vars],0,'post'), -ones(num_vars,1)]};
            end
            return
        end
        
        %b is not estimated
        if num_ud == 0
            A_Bu = mrdivide(x_seg_plus,[x_seg]);
            dA = A_Bu(:,1:num_vars);
            dB = zeros(size(A_Bu,1),1);
            dC = eye(num_vars); dD = zeros(num_vars, num_ud+1);
            sysd = ss(dA, dB, dC, dD, Ts);
            try
                sysc = d2c(sysd,"zoh");
            catch
                continue; %try with less derivatives included
            end
            ode(label) = {[padarray(sysc.A,[max_num_var-num_vars max_num_var-num_vars],0,'post'), padarray(sysc.B,[max_num_var-num_vars],0,'post')]};
            break; %ode found
        %b is estimated
        else
            A_Bu = mrdivide(x_seg_plus,[x_seg;ud_seg]);
            dA = A_Bu(:,1:num_vars);
            dB = A_Bu(:,num_vars+1:end);
            dC = eye(num_vars); dD = zeros(num_vars, num_ud+1);
            sysd = ss(dA, dB, dC, dD, Ts);
            try
                sysc = d2c(sysd,"zoh");
            catch
                continue; %try with a less derivatives included
            end
            ode(label) = {[padarray(sysc.A,[max_num_var-num_vars max_num_var-num_vars],0,'post'), padarray(sysc.B,[max_num_var-num_vars],0,'post')]};
            break; %ode found
        end
    end
end