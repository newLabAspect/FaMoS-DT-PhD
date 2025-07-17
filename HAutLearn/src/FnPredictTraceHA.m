function [sim_trace] = FnPredictTraceHA(trace,conditions,ode)
% FnPredictTraceHA uses the HA-model of the system and inital conditions of
% a preexisting trace to create a trace prediction
    global num_var num_ud tolLI Ts offsetCluster

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
        state_id = find(trace.labels_num == curr_state);
        A = cell2mat(ode(state_id));
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
        % Check if a transitions condition holds, thus transition needs to be taken
        for j = 1:size(conditions,1)
            condition = conditions(j,:);
            % Considered condition irrelevant if origin is not current state
            if(condition(1,1) == curr_state)
                % Evaluate considered condition
                overall_sum = +1.0;
                for k = 1:num_var
                    overall_sum = overall_sum + new_x(k,1) * condition(1,2+k);
                end
                for k = 1:num_ud
                    overall_sum = overall_sum + sim_ud(i,k) * condition(1,2+num_var+k);
                end
                % Considered condition holds, thus switch states
                if((overall_sum < tolLI && condition(1,2+num_var+num_ud+1) == -1) || ...
                   (overall_sum > -tolLI && condition(1,2+num_var+num_ud+1) == +1))
                    curr_state = condition(1,2);
                    break;
                end
            end
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