function [correct, false] = FnEvalCluster(estimated, truth, truth_chps)
% FnEvalCluster evalautes the generated clusters against the ground truth 
% clusters provided in truth. Because truth includes segments which are 
% too short to detect, all changepoints of truth are provided to filter these  
    global windowSize

    % Remove too small segments out of truth array
    indx_remove = find((truth_chps(2:end) - truth_chps(1:(length(truth_chps)-1))) <= windowSize);
    truth(indx_remove) = [];

    % Find fitting permutation between generated and truth cluster ids.
    % Generate mapping truth to generated cluster ids, but this masks if  
    % there are not enough generated clusters  
    labels_truth = unique(truth);
    M_T2E = containers.Map('KeyType','double','ValueType','double');
    for i = 1:length(labels_truth)
        % Find all segments in truth with i'th cluster id
        indx = find(truth==labels_truth(i));
        % Remove indexes that are out of bounds for generated
        indxOutOfBounds = find(indx > size(estimated,1));
        indx(indxOutOfBounds) = [];
        % Associate generated cluster id which occurs most often at these segments
        M_T2E(labels_truth(i)) = mode(estimated(indx));
    end

    % Generate mapping generated to truth cluster ids, but this masks if  
    % there are too many generated clusters  
    labels_est = unique(estimated);
    M_E2T = containers.Map('KeyType','double','ValueType','double');
    for i = 1:length(labels_est)
        % Find all segments in generated with i'th cluster id
        indx = find(estimated==labels_est(i));
        % Remove indexes that are out of bounds for truth
        indxOutOfBounds = find(indx > size(truth,1));
        indx(indxOutOfBounds) = [];
        % Associate truth cluster id which occurs most often at these segments
        M_E2T(labels_est(i)) = mode(truth(indx));
    end

    % Evaluate sequence of generated cluster ids against ground truth using
    % both mappings to overcome their respective limitations
    correct = 0;
    false = 0;
    for i = 1:max(size(estimated,1),size(truth,1))
        if(i > size(truth,1) || i > size(estimated,1))
            % Too many/few segments in sequence of generated cluster ids
            % Exclude trailing zeros as they are placeholders
            if(i <= size(estimated,1) && estimated(i,1) ~= 0)
                false = false + 1;
            end
        elseif(estimated(i,1) == M_T2E(truth(i,1)) && M_E2T(estimated(i,1)) == truth(i,1))
            % Both mappings indicate correct result
            correct = correct + 1;
        else
            false = false + 1;
        end
    end
end