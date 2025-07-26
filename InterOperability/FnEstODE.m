function [ode, all_labels] = FnEstODE(trace)
% FnEstODE estimates the ODE for each cluster in the trace datastructure
global  num_var num_ud offsetCluster

% Iterate over all clusters
all_labels = [];
for j = 1:length(trace)
    all_labels = union(all_labels,unique(trace(j).labels_trace));
end
all_labels = sort(all_labels);
len_labels = length(all_labels);
num_vars = num_var*(1+offsetCluster); 
for label = 1:len_labels
    x_seg = [];
    x_seg_plus = [];
    ud_seg = [];
    % Extract all segments from all traces associated with current
    % cluster id to save state (and input) vector
    for j = 1:length(trace)
        idx = find(trace(j).labels_trace == all_labels(label));
        x = trace(j).xs(:,1:num_vars);
        ud = trace(j).ud(:,1:num_ud);
         
        startj = trace(j).chpoints(idx);
        endj = trace(j).chpoints(idx+1)-1;
        
        for n = 1:length(startj)
            % Save state space vector
            x_seg = [x_seg, x(startj(n):(endj(n)-1), :)'];
            % Additional input with only 1s considered to allow for
            % additive constants
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

    % ODE datastructure is A appended by B. As different ODE degrees
    % are supported, A and B are extended to their maximal sizes
    
    % No conversion to a continuous model, denoted by -1s in ODE
    A_Bu = mrdivide(x_seg_plus,[x_seg;ud_seg]);
    %dA = A_Bu(:,1:num_vars);
    %dB = A_Bu(:,num_vars+1:end);
    ode(label) = {A_Bu}; %{[dA, dB]};
end
end