function [mse, sim_trace] = FnEvalFlowAccuracy(trace, ode)
%FnEvalFlowAccuracy Determines the accuracy of the flow function on a
%trace.
%   Uses the identified segments and modes to evaluate the acuracy of the
%   flow functions
     global useTime num_var num_ud Ts offsetCluster

    % Preallocated arrays for speed
    sim_x = zeros(size(trace.x,1),num_var*(1+offsetCluster));
    sim_state = zeros(size(trace.x,1),1);
    sim_chp = [];
    sim_labels = [];

    % Extract initial conditions at predefined position from trace datastructure
    offsetPred = 0+1;
    sim_x(1,:) = trace.xs(offsetPred,1:num_var*(1+offsetCluster));
    indxLastSwitch = find(trace.chpoints(:,1) <= offsetPred,1,'last');
    curr_state = trace.labels_trace(indxLastSwitch,1);
    % Inputs are fully needed for predicition
    sim_ud = trace.ud;
    chpoints = trace.chpoints;
    mse = zeros(num_var*(1+offsetCluster),1);

    % Trace Prediction (of timepoint i using timepoint i-1)
    for j = 1:length(chpoints)-1
        for i = chpoints(j):(chpoints(j+1)-1)
            % Predict next data point based on ODE
            A = cell2mat(ode(curr_state));
            B = A(1:size(A,1),(size(A,1)+1):end);
            A = A(1:size(A,1),1:size(A,1));
            curr_x = trace.xs(i,1:num_var*(1+offsetCluster))';
            curr_u = [1];
            if num_ud ~= 0
                curr_u = [trace.ud(i,1:num_var*(1+offsetCluster))';1];
            end
            new_x = A * curr_x + B * curr_u;
            sim_x(i,:) = new_x';
        end
    end

    for j = 1:length(mse)
        mse(j) = mean((sim_x(:,j)-trace.xs(:,j)).^2);
    end

    % Create trace data structure
    sim_trace.x = sim_x;
    sim_trace.ud = sim_ud;
end