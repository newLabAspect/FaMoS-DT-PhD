function [sim_trace] = FnPredictTraceDTL(trace,Mdl,ode)
% FnPredictTraceDTL uses the DT-model of the system and inital conditions of
% a preexisting trace to create a trace prediction
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
    lastSwitch = trace.chpoints(indxLastSwitch,1);
    sim_state(1,1) = trace.labels_trace(indxLastSwitch,1);
    curr_state = trace.labels_trace(indxLastSwitch,1);
    % Inputs are fully needed for predicition
    sim_ud = trace.ud;

    % Trace Prediction (of timepoint i using timepoint i-1)
    for i = (offsetPred+1):size(trace.x)
        % Predict next data point based on ODE
        A = cell2mat(ode(curr_state));
        B = A(1:size(A,1),(size(A,1)+1):end);
        A = A(1:size(A,1),1:size(A,1));
        curr_x = sim_x(i-1,1:size(A,1))';
        curr_u = [1];
        if num_ud ~= 0
            curr_u = [sim_ud(i,:)';1];
        end
        new_x = A * curr_x + B * curr_u;
        sim_x(i,:) = new_x';

        % Predict next system mode
        sim_state(i,1) = curr_state;
        last_state = curr_state;
        % Construct feature vector for timepoint i-1
        if (num_ud ~= 0)
            values = [curr_x(1:num_var,1)' sim_ud(i,1:num_ud)];
        else
            values = [curr_x(1:num_var,1)'];
        end
        if (useTime)
            curr_state = predict(Mdl,[last_state values (i-1)-lastSwitch]); 
        else
            curr_state = predict(Mdl,[last_state values]);
        end

        % Track if state has changed (needed to compute time since last switch)
        if(last_state ~= curr_state)
            lastSwitch = i;
            sim_chp = [sim_chp; i]; % maybe i+1 or i-1 
            sim_labels = [sim_labels; last_state];
        end
    end

    % Create trace data structure
    sim_trace.x = sim_x;
    sim_trace.ud = sim_ud;
    sim_trace.chpoints = [1; sim_chp; length(sim_x)];
    sim_trace.labels_trace = sim_labels;
    sim_trace.labels_num = unique(sim_labels);
end