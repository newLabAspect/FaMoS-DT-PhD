function [trace] = FnCleanChangePoints(trace)
%FNCLEANCHANGEPOINTS Summary of this function goes here
%   Detailed explanation goes here
    
    for i = 1:length(trace)
        indx = find(trace(i).labels_trace(1:(end-1),1) == trace(i).labels_trace(2:end,1));
        if isempty(indx)
            continue;
        end
        trace(i).labels_trace(indx+1) = [];
        del_chps = trace(i).chpoints(indx+1);
        trace(i).chpoints(indx+1) = [];
        for j = 1:length(trace(i).chpoints_per_var)
            chps = cell2mat(trace(i).chpoints_per_var(j,1));
            [~,indx] = ismember(del_chps,chps);
            indx = find(indx>0);
            if isempty(indx)
                continue;
            end
            chps(indx) = [];
            trace(i).chpoints_per_var(j,1) = {chps};
        end
    end
end

