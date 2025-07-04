function [omega_seg] = FnEvalChangePoints(chps, truth_chps)
%FNEVALCHANGEPOINTS Calculate deviation between predicted and ground truth switches
%
% Input:
%   chps - array of predicted switch points
%   truth_chps - array of ground truth switch points
%
% Output:
%   deviation - calculated deviation score
    
    length_penalty = 10;
    
    deviation = 0.0;
    
    % Add length penalty
    deviation = deviation + length_penalty * abs(length(truth_chps) - length(chps));
    
    % Initialize matching array
    matching = nan(length(truth_chps), 1);
    
    % Find best matches for each ground truth switch
    for i = 1:length(truth_chps)
        % Find the switch with minimum distance
        distances = abs(chps - truth_chps(i));
        [~, min_idx] = min(distances);
        matching(i) = min_idx;
        
        % Check if this index was already used
        already_used = false;
        used_idx = 0;
        for j = 1:(i-1)
            if ~isnan(matching(j)) && matching(j) == min_idx
                already_used = true;
                used_idx = j;
                break;
            end
        end
        
        if already_used
            % Calculate distances for conflict resolution
            new_dist = abs(chps(min_idx) - truth_chps(i));
            old_dist = abs(chps(min_idx) - truth_chps(used_idx));
            
            if new_dist > old_dist
                matching(i) = NaN;
            else
                matching(used_idx) = NaN;
            end
        end
    end
    
    % Calculate final deviation
    for i = 1:length(matching)
        if ~isnan(matching(i))
            deviation = (deviation + abs(chps(matching(i)) - truth_chps(i))) / length(chps);
        end
    end
    omega_seg = deviation;
end
