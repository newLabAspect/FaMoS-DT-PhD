function [correct, false] = FnEvalCluster(estimated, truth, truth_chps)
    global windowSize
    %remove too small segments out of truth array; TODO: Which threshold from codebase to use
    indx_remove = find((truth_chps(2:end) - truth_chps(1:(length(truth_chps)-1))) <= windowSize);
    truth(indx_remove) = [];

    %find fitting permutation, both directions used because, e.g., T2E masks 
    %if there are not enough states in estimated (and vice-versa)  
    labels_truth = unique(truth);
    M_T2E = containers.Map('KeyType','double','ValueType','double');
    for i = 1:length(labels_truth)
        indx = find(truth==labels_truth(i));
        %remove indx that are out of bounds for truth array
        indxOutOfBounds = find(indx > size(estimated,1));
        indx(indxOutOfBounds) = [];
        M_T2E(labels_truth(i)) = mode(estimated(indx));
    end

    labels_est = unique(estimated);
    M_E2T = containers.Map('KeyType','double','ValueType','double');
    for i = 1:length(labels_est)
        indx = find(estimated==labels_est(i));
        %remove indx that are out of bounds for truth array
        indxOutOfBounds = find(indx > size(truth,1));
        indx(indxOutOfBounds) = [];
        M_E2T(labels_est(i)) = mode(truth(indx));
    end

    %decide correct false incl. some edge cases
    correct = 0;
    false = 0;
    for i = 1:max(size(estimated,1),size(truth,1))
        if(i > size(truth,1) || i > size(estimated,1))
            %sometimes estimated filled with zeros at end
            if(i <= size(estimated,1) && estimated(i,1) ~= 0)
                false = false + 1;
            end
        elseif(estimated(i,1) == M_T2E(truth(i,1)) && M_E2T(estimated(i,1)) == truth(i,1))
            correct = correct + 1;
        else
            false = false + 1;
        end
    end
end