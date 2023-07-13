function [sim_x, sim_state] = FnPredictTraceDTL(trace,Mdl,ode)
% FnPredictTraceDTL uses the DT-model of the system and inital conditions of
% a preexisting trace to create a trace prediction
    global useTime num_var Ts precisionDTL

    % Preallocated arrays for speed
    sim_x = zeros(size(trace.x,1),size(trace.x,2));
    sim_state = zeros(size(trace.x,1),1);

    % Extract initial conditions at predefined position from trace datastructure
    offsetPred = 0+1;
    sim_x(1,:) = trace.x(offsetPred,:);
    indxLastSwitch = find(trace.chpoints(:,1) <= offsetPred,1,'last');
    lastSwitch = trace.chpoints(indxLastSwitch,1);
    sim_state(1,1) = trace.labels_trace(indxLastSwitch,1);
    curr_state = trace.labels_trace(indxLastSwitch,1);
    
    % Trace Prediction (of timepoint i using timepoint i-1)
    for i = (offsetPred+1):size(trace.x)
        % Predict next data point based on ODE
        A = cell2mat(ode(curr_state));
        A = A(1:size(A,1),1:size(A,1));
        curr_x = sim_x(i-1,1:size(A,1))';
        new_x_dot = A * curr_x;
        % Assume: First rows of A represent integration, thus leave out
        new_x = [curr_x(1:end,1); new_x_dot((end-num_var+1):end,1)];
        % Update all derivatives bottom-up with the rectangle integration rule
        for j = (size(new_x,1)-num_var):-1:1
            new_x(j,1) = new_x(j,1) + new_x(j+num_var,1) * Ts;
        end
        % Pad with zeros as constant row length needed
        sim_x(i,:) = [new_x', zeros(1,size(sim_x,2)-size(new_x,1))];

        % Predict next system mode
        sim_state(i,1) = curr_state;
        last_state = curr_state;
        % Construct feature vector for timepoint i-1
        if (useTime)
            curr_state = predict(Mdl,[last_state FnRoundToInterval(curr_x(1:num_var,1)',precisionDTL) (i-1)-lastSwitch]); 
        else
            curr_state = predict(Mdl,[last_state FnRoundToInterval(curr_x(1:num_var,1)',precisionDTL)]);
        end

        % Track if state has changed (needed to compute time since last switch)
        if(last_state ~= curr_state)
            lastSwitch = i;
        end
    end
end