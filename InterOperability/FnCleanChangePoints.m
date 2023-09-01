function [trace] = FnCleanChangePoints(trace)
% FnCleanChangePoints removes unnecessary changepoints between segments
% with the same cluster id
%   These are mainly caused due to changes in the input traces that affect
%   the dynamic behavior of the system but does not trigger a mode switch
    
    for i = 1:length(trace)
        % Find adjacent segments with the same cluster id 
        indx = find(trace(i).labels_trace(1:(end-1),1) == trace(i).labels_trace(2:end,1));
        if isempty(indx)
            continue;
        end
        % Remove intermediate changepoint (globally) and repition of label
        trace(i).labels_trace(indx+1) = [];
        del_chps = trace(i).chpoints(indx+1);
        trace(i).chpoints(indx+1) = [];
        % Remove intermediate changepoint (locally) on each variable
        for j = 1:length(trace(i).chpoints_per_var)
            chps = cell2mat(trace(i).chpoints_per_var(j,1));
            [~,indx] = ismember(del_chps,chps);
            % Zero indicates that no occurence was found, thus remove from
            % index vector
            indx = find(indx>0);
            if isempty(indx)
                continue;
            end
            chps(indx) = [];
            trace(i).chpoints_per_var(j,1) = {chps};
        end
    end
end

