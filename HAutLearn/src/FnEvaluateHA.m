function [correct,false] = FnEvaluateHA(trace, conditions, tol)
% FnEvaluateHA evaluates the trained HA model against an evaluation dataset
% provided by a trace datastructure
global num_var
    
    correct = 0;
    false = 0;
    % Mutliple traces can be used in the evaluation
    for l = 1:length(trace)
        state_index = 1;
        % All datapoints with an successor are evaluated
        for i = 1:(size(trace(l).x,1)-1)
            % Changepoint is reached, thus next segment has to be considered
            if(i>trace(l).chpoints(state_index+1)) 
                state_index = state_index + 1;
            end

            % Compute current and next predicted state
            state = trace(l).labels_trace(state_index);
            predicted_state = state;
            % Check if a transitions condition holds, thus transition needs to be taken
            for j = 1:size(conditions,1)
                condition = conditions(j,:);
                % Considered condition irrelevant if origin is not current state
                if(condition(1,1) == state)
                    % Evaluate considered condition
                    overall_sum = +1.0;
                    for k = 1:num_var
                        overall_sum = overall_sum + trace(l).x(i,k) * condition(1,2+k);
                    end
                    % Considered condition holds, thus switch states
                    if((overall_sum < tol && condition(1,2+num_var+1) == -1) || ...
                       (overall_sum > -tol && condition(1,2+num_var+1) == +1))
                        predicted_state = condition(1,2);
                        break;
                    end
                end
            end

            % Compare predicted next state to actual next state
            next_state = state;
            if(i+1>trace(l).chpoints(state_index+1))
                next_state = trace(l).labels_trace(state_index+1);
            end
            % Update respective counter
            if(next_state == predicted_state)
                correct = correct + 1;
            else
                false = false + 1;
                %disp([num2str(overall_sum),int2str(condition(1,2+num_var+1))])
            end
        end
    end
end