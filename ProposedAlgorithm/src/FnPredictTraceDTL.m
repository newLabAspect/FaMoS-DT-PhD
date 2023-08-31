function [sim_trace] = FnPredictTraceDTL(trace,Mdl,ode)
% FnPredictTraceDTL uses the DT-model of the system and inital conditions of
% a preexisting trace to create a trace prediction
    global useTime num_var num_ud Ts

    % Preallocated arrays for speed
    sim_x = zeros(size(trace.x,1),size(trace.x,2));
    sim_state = zeros(size(trace.x,1),1);
    sim_chp = [];
    sim_labels = [];

    % Extract initial conditions at predefined position from trace datastructure
    offsetPred = 0+1;
    sim_x(1,:) = trace.x(offsetPred,:);
    sim_ud = trace.ud;
    indxLastSwitch = find(trace.chpoints(:,1) <= offsetPred,1,'last');
    lastSwitch = trace.chpoints(indxLastSwitch,1);
    sim_state(1,1) = trace.labels_trace(indxLastSwitch,1);
    curr_state = trace.labels_trace(indxLastSwitch,1);

    % Remove zero rows/columns in ODEs (needed for earlier HAutLearn compatibility)
%     for k = 1:length(ode)
%         A = cell2mat(ode(k));
%         A = shrinkMatrix(A);
%         ode(k) = {A};
%     end
    
    % Trace Prediction (of timepoint i using timepoint i-1)
    for i = (offsetPred+1):size(trace.x)
        % Predict next data point based on ODE
        A = cell2mat(ode(curr_state));
        B = A(1:size(A,1),(size(A,1)+1):end);
        A = A(1:size(A,1),1:size(A,1));
        curr_x = sim_x(i-1,1:size(A,1))';
        new_x_dot = A * curr_x + B * [sim_ud(i,:)';1];

        % Assume: First rows of A represent integration, thus leave out
        new_x = [curr_x(1:end,1); new_x_dot((end-num_var+1):end,1)];
        % Update all derivatives bottom-up with the rectangle integration rule
        for j = (size(new_x,1)-num_var):-1:1
            new_x(j,1) = new_x(j,1) + new_x(j+num_var,1) * Ts;
        end
        % Pad with derivatives as constant row length needed
        sim_x(i,:) = [new_x', zeros(1,size(sim_x,2)-size(new_x,1))];
        for j = (size(new_x,1)+1):size(sim_x,2)
            sim_x(i,j) = (sim_x(i,j-num_var) - sim_x(i-1,j-num_var)) / Ts;
        end

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

function A_new = shrinkMatrix(A)
% sprinkMatrix reduces the matrix A to the minimal quadratic size by
% deleting zero columns and rows
    
    % Find all non-zero entries
    [row,col] = find(A);
    % Determine minimal quadratic size
    len = max(max(row),max(col));
    % Return result
    A_new = A(1:len,1:len);
end