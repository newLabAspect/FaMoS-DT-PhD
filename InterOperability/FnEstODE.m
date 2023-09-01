function ode = FnEstODE(trace)
% FnEstODE estimates the ODE for each cluster in the trace datastructure
global Ts num_var num_ud offsetCluster

% Iterate over all clusters
len_labels = length(trace(1).labels_num);
for label = 1:len_labels
    max_num_var = num_var * (1 + offsetCluster);
    % If higher order ODEs are estimated, derivatives need to be included
    % in the state vector. But ODE degree unknown, thus searched greedy
    for off = offsetCluster:-1:0
        x_seg = [];
        x_seg_plus = [];
        ud_seg = [];
        num_vars = num_var * (1+off);
        % Extract all segments from all traces associated with current
        % cluster id to save state (and input) vector
        for j = 1:length(trace)
            labels_trace = trace(j).labels_trace;
            idx = find(labels_trace == label);
            x = trace(j).x(:,1:num_vars);
            ud = trace(j).ud;
             
            startj = trace(j).chpoints(idx);
            endj = trace(j).chpoints(idx+1)-1;
            
            for n = 1:length(startj)
                % Save state space vector
                x_seg = [x_seg, x(startj(n):(endj(n)-1), :)'];
                % If no input available, use 1s as placeholders
                if isempty(ud)
                    ud_seg = [ud_seg, ones(1, (endj(n)- startj(n)))];
                else
                    ud_seg = [ud_seg, [ud(startj(n):(endj(n)-1), :)'; ones(1, (endj(n)- startj(n)))]];
                end
                % Save state space vector delayed by one timestep as 
                % discrete state space representation used
                x_seg_plus = [x_seg_plus, x((startj(n)+1):endj(n), :)'];
            end
        end

        % ODE datastructure is A appended by b. As different ODE degrees
        % are supported, A and b are extended to their maximal sizes
        
        % No conversion to a continous model possible, denoted by -1s in ODE
        if size(x_seg,2)/length(trace)==2
            A_Bu = mrdivide(x_seg_plus,[x_seg;ud_seg]);
            dA = A_Bu(:,1:num_vars);
            dB = A_Bu(:,num_vars+1:end);
            % -1 indicate the ode is descrete and is used for variable reset
            ode(label) = {[fillToSize(dA,[max_num_var max_num_var]), fillToSize(dB,[max_num_var 1]), -ones(num_vars,1)]};
            break;
        end
        
        A_Bu = mrdivide(x_seg_plus,[x_seg;ud_seg]);
        dA = A_Bu(:,1:num_vars);
        dB = A_Bu(:,num_vars+1:end);
        dC = eye(num_vars); dD = zeros(num_vars, num_ud+1);
        sysd = ss(dA, dB, dC, dD, Ts);
        try
            sysc = d2c(sysd,"zoh");
        catch
            % Conversion to continous model not possible because too
            % many derivatives included in state vector, retry with less
            continue;
        end
        ode(label) = {[fillToSize(sysc.A,[max_num_var max_num_var]), fillToSize(sysc.B,[max_num_var 1])]};
        break;
    end
end
end

function A_new = fillToSize(A,len)
% fillToSize extends the matrix A to size len x len
    A_new = zeros(len);
    A_new(1:size(A,1),1:size(A,2)) = A;
end