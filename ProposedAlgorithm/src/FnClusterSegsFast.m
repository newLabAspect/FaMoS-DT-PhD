function trace = FnClusterSegsFast(trace, x, ud)
% FnClusterSegsFast computes the clusters of the trace using DTW-comparisons 
% and uses LMI-comparisons to refine the results if selected

    % Compute local and global segments from chpoints in trace datastructure
    [segIndex, segIndex_var] = computeSegments(trace);
    
    % Compute similarity matrix which is used for clustering in next step
    combined_metric = computeSimilarityMatrix(x,segIndex_var);

    % Create local clusters on each output variable
    [cluster_segs, trace] = computeClustersLocal(trace, combined_metric, segIndex_var);

    % Merge local clusters to global clusters, potentially refine with LMI
    [cluster_global, trace] = computeClustersGlobal(x, trace, cluster_segs, segIndex, segIndex_var);

    % Save global results into trace data structure
    labels_num = unique(cluster_global(:,1));

    for i=1:length(trace)
        chpoints = (trace(i).chpoints);
        len_segs = length(chpoints)-1;
        trace(i).labels_num = labels_num;
        trace(i).labels_trace = cluster_global(1:len_segs,1);
        cluster_global(1:len_segs,:) = [];
    end
end

function [segIndex, segIndex_var] = computeSegments(trace)
% computeSegments returns local (per output variable) and global segments
% based on local and global changepoints saved in trace data structure
    global num_var
    % Include dummy entry to simplify computations
    segIndex = [0,0];
    segIndex_var = cell(num_var,1);
    segIndex_var(1:num_var,1) = {[0,0]};
    % Transform changepoints to segments consisting of start and end point
    for i=1:length(trace)  
        % Global segments
        chpoints = (trace(i).chpoints);
        chsegments = [chpoints(1:end-1), chpoints(2:end)];
        % Increase start point by 1 to have no overlapping segments
        chsegments(2:end,1) = chsegments(2:end,1)+1;
        % Apply offsets to changepoints as all trace values are appended to 
        % form x (and ud) to ensure consistency
        segIndex = [segIndex; segIndex(end,2) + chsegments];

        % Local segments (per output variable)
        chp_var = (trace(i).chpoints_per_var);
        for j=1:num_var
            chp_curr = cell2mat(chp_var(j,1));
            chsegments_var = [chp_curr(1:end-1),chp_curr(2:end)];
            % Increase start by 1 as done above
            chsegments_var(2:end,1) = chsegments_var(2:end,1)+1;
            % Apply offsets to be consistent with x (and ud) as above
            segIndex_var_temp = cell2mat(segIndex_var(j,1));
            segIndex_var(j,1) = {[segIndex_var_temp; segIndex_var_temp(end,2) + chsegments_var]};
        end 
    end
    % Dummy lines are deleted (used so that first append offset is 0)
    segIndex(1,:) = [];
    for j=1:num_var
        segIndex_var_temp = cell2mat(segIndex_var(j,1));
        segIndex_var_temp(1,:) = [];
        segIndex_var(j,1) = {segIndex_var_temp};
    end
end

%% Functions related to Computation of the similarity matrix

function sim = computeComparison(seg_1,seg_2)
% computeComparison returns the similarity index of two segments
%   The similarity index is computed using a DTW comparison and combining
%   the distance and diagonality metric into a single one
    [dist, i_x, i_y] = dtw(seg_1,seg_2);
    diag = corrcoef(i_x,i_y);
    diag = diag(1,2);
    dist = dist/size(i_x,2);
    % Diag can be NaN if one signal is a constant resulting in
    % zero variance, in this case no similarity is given
    if isnan(diag)
        diag = 0.0;
    end
    % Combine both metrics into a single one
    sim = 0.5*dist+0.5*(1-diag);
end

function combined_metric = computeSimilarityMatrix(x, segIndex_var)
% computeSimilarityMatrix returns the similarity matrix which quantifies
% how similar each two segments are with repect to each output variable for 
% the selected derivative
%   The ideas to shorten the segments to a common length and carry out two
%   comparisons - one aligned at the end and one aligned at the start - are
%   used to compute the similarity matrix
    global num_var winlen offsetCluster

    combined_metric = cell(num_var,1);
    % Comparison Metric computed for each output variable
    for k = 1:num_var
        segIndex_curr = cell2mat(segIndex_var(k,1));
        combined_curr = zeros(size(segIndex_curr,1),size(segIndex_curr,1));
        % Pairwise comparison first index
        for i = 1:size(segIndex_curr,1)
            % Start and end trimmed to ignore, e.g., peaks in derivatives
            % caused by changepoints
            start_i = segIndex_curr(i,1)+winlen;
            end_i = segIndex_curr(i,2)-winlen;
            seg_i = x(start_i:end_i,k+offsetCluster);
            % Pairwise comparison second index
            for j = i:size(segIndex_curr,1)
                start_j = segIndex_curr(j,1)+winlen;
                end_j = segIndex_curr(j,2)-winlen;
                seg_j = x(start_j:end_j,k+offsetCluster);
                common_len = min(size(seg_i,1),size(seg_j,1));
                
                % Constant inital offset removed if original trace value used
                incConst = 1;
                if offsetCluster ~= 0
                    incConst = 0;
                end
                
                % Compute DTW-Comparison for shortened, at END aligned signals
                sim_index_end = computeComparison(seg_i((end-common_len+1):end,1)-seg_i(end,1)*incConst,seg_j((end-common_len+1):end,1)-seg_j(end,1)*incConst);
                
                % Compute DTW-Comparison for shortened, at START aligned signals
                sim_index_start = computeComparison(seg_i(1:common_len,1)-seg_i(1,1)*incConst,seg_j(1:common_len,1)-seg_j(1,1)*incConst);

                % Select the better similarity index and save it in
                % symmetrical extended matrix
                if (sim_index_start <= sim_index_end)
                    combined_curr(i,j) = sim_index_start;
                    combined_curr(j,i) = sim_index_start;
                else
                    combined_curr(i,j) = sim_index_end;
                    combined_curr(j,i) = sim_index_end;
                end
            end
        end
        % Add matrix of current variable to general matrix
        combined_metric(k,1) = {combined_curr};
    end

end

%% Functions related to clustering

function [clusters,alreadyClustered] = FnDecideSimilar(i,j,clusters,alreadyClustered, combined_metric, thres)
% FnDecideSimilar decides if two segments are similar and updates the
% clusters if this is the case
    similar = true;
    
    % Segments are not similar, thus they should not be in the same cluster
    if(combined_metric(i,j) > thres)
        similar = false;
    end

    % Segments are similar and they are similar at highest confidence level
    % analyzed up till now (thus lowest similarity index up till now)
    if(similar && (alreadyClustered(j,1) == 0 || alreadyClustered(j,2) > combined_metric(i,j)))
        % Second segment was clustered before, but at a lower confidence level
        % => remove it from previous cluster (these clusters will be merged over time)
        if(alreadyClustered(j,1) ~= 0)
            clusters(alreadyClustered(j,1)) = {setdiff(cell2mat(clusters(alreadyClustered(j,1))),j)};
        end
        % First segment was clustered before => add second segment to the 
        % cluster of the first segment
        if(alreadyClustered(i,1) ~= 0)
            clusters(alreadyClustered(i,1)) = {union(cell2mat(clusters(alreadyClustered(i,1))),j)};
            alreadyClustered(j,1) = alreadyClustered(i,1);
        % First segment was not clustered before => create new cluster with id i
        else
            clusters(i) = {union(cell2mat(clusters(i)),j)};
            alreadyClustered(j,1) = i;
        end
        % Save confidence level for later comparisons
        if(i ~= j)
            alreadyClustered(j,2) = combined_metric(i,j);
        end
    end
end

% Carried over from HAutLearn with minor modifications as this function
% should only provide the result similar or not
function result = FnDecideSimilarLMI(xseg1,udseg1,xsegj,udsegj)  
    global sigma winlen
       
    %variable number
    num_var = size(xseg1,1);
    %size of ud
    num_ud = size(udseg1,1);
    
    %xk1: delete some points in end bc of inaccuracy of chpoints
    xk1 = xseg1(:,1:end - winlen); 
    ud1 = [];
    if num_ud ==0
        %not used: ud1 has same dimensions as xk1 but all are ones
        ud1 = ones(1,size(xk1,2));
    else
        ud1 = udseg1(:,1:end - winlen);
        ud1 = [ud1; ones(1,size(xk1,2))];
    end
    %xki1: delete some points in beginning bc of chpoints inaccuracy
    xki1 =  xseg1(:,winlen+1:end);
    
    %Ok is xk1 in the first row and ud1 in the 2nd row (prep LMI)
    Ok = [xk1;ud1];
    Oki = xki1;
    
    setlmis([]); % Set a new LMI
    structV = [1, 1+num_ud+num_var]; %matrix structure
    for i = 1:num_var % create Vi matrix for each variable
        V{i} = lmivar(2, structV);
    end
    
    %setup up Linear Matrix Inequalities
    error_t = (sigma * size(Oki,2))^2; %error tolerance
    for i = 1:size(Oki,1)
        n = (i+(1-1)*num_var);
        lmiterm([-n, 1,1,0],1);
        lmiterm([-n, 2,1,V{i}],1,Ok)
        lmiterm([-n, 2,1,0],-Oki(i,:))
        lmiterm([-n, 2,2,0], error_t)
    end
    lmi = getlmis; 
    
    %xkj eq. to xk1, udj eq. ud1, xkij eq. xki1, Okj eq. Ok and Okij eq. Oki
        
    xkj = xsegj(:,1:end - winlen);
    xkij = [];
    xkij = xsegj(:, winlen+1:end);
    
    if num_ud == 0
        udj = ones(1,size(xkj,2));
    else
        udj = udsegj(:,1:end - winlen);
        udj = [udj; ones(1,size(xkj,2))];
    end

    Okj = [xkj; udj];
    Okij = xkij;
    
    %append to already existing LMI
    setlmis(lmi);
    error_tj = (sigma * size(Okj,2))^2;
    for i = 1:num_var
        n = (i+(2-1)*num_var);
        lmiterm([-n, 1,1,0],1);
        lmiterm([-n, 2,1,V{i}],1,Okj)
        lmiterm([-n, 2,1,0],-Okij(i,:))
        lmiterm([-n, 2,2,0], error_tj)
    end
    lmi = getlmis; 
    newlmi = lmi;  
    %solve LMI
    [tmin,xfeas] = feasp(newlmi);
    %v = dec2mat(Ilmi,xfeas,V);
    %1 and j are similar
    if (tmin <=0)     
        result = true;
        return
    %1 and j are not similar
    else
        result = false;
        return
    end
end

function cluster_global = refineClustersLMI(x,cluster_global,segIndex)
% refineClustersLMI checks all cluster pairs for similarity using an
% LMI-based approach and merges them if they are similar
%   Currently the segment which occurs first and is associated with that
%   cluster id is used for the LMI comparison
    global num_var winlen offsetCluster windowSize

    labels_num = unique(cluster_global(:,1));
    progress_cnt = 1;

    % Check for each cluster pair if they are similar, first entry
    for i=labels_num'
        % Select segment which is part of first cluster (more advanced selection?)
        posi = find(cluster_global == i,1);
        % Check if cluster still existent or if it was already merged
        if(isempty(posi))
            continue;
        end
        start_i = segIndex(posi,1)+winlen;
        end_i = segIndex(posi,2)-winlen;
        % Global clustering only possible if segment long enough (tune factor?)
        if(end_i-start_i-1 < 3 * windowSize)
            continue;
        end
        % Consider derivatives up to selected degree as done in pure LMI approach
        seg_i = x(start_i:end_i,1:(num_var+offsetCluster*num_var));

        % Check for each cluster pair if they are similar, second entry
        for j=setdiff(labels_num',labels_num(1:i)')
            disp([int2str(progress_cnt),'/', int2str((length(labels_num)^2-length(labels_num))/2)])
            progress_cnt = progress_cnt + 1;
            % Extract seg_j analog to seg_i before
            posj = find(cluster_global == j,1);
            if(isempty(posj))
                continue;
            end
            start_j = segIndex(posj,1)+winlen;
            end_j = segIndex(posj,2)-winlen;
            if(end_j-start_j-1 < 3 * windowSize)
                continue;
            end
            seg_j = x(start_j:end_j,1:(num_var*(offsetCluster+1)));

            % If clusters are similar, merge them
            res = FnDecideSimilarLMI(seg_i',[],seg_j',[]);
            if(res)
                posk = find(cluster_global == j);
                cluster_global(posk,1) = i;
            end
        end
    end
end

function [cluster_segs, trace] = computeClustersLocal(trace, combined_metric, segIndex_var)
% computeClustersLocal returns clusters computed on each output variable using the 
% similarity matrix with an threshold-based approach with dynamically computed thresholds
    global num_var thresClusterMin thresClusterMax facThres
    % Extend trace data structure to save sequence of clusters for each var 
    for i=1:length(trace)
        trace(i).labels_trace_per_var = cell(num_var,1);
    end

    % Carry out clustering on each output variable using the corresponding
    % set of local changepoints
    cluster_segs = cell(num_var,1);
    for k = 1:num_var
        segIndex_curr = cell2mat(segIndex_var(k,1));
        clusters = cell(size(segIndex_curr,1),1);
        alreadyClustered = zeros(size(segIndex_curr,1),2);
        for i = 1:size(segIndex_curr,1)
            % Due to assumed sequential behavior of HA, at least one other
            % segment is similar, thus threshold computed based on it
            thres = facThres*min(combined_metric{k,1}(i,union(1:(i-1),(i+1):end)));
            % But threshold limited to an interval
            if(thres < thresClusterMin)
                thres = thresClusterMin;
            elseif(thres > thresClusterMax)
                thres = thresClusterMax;
            end
            % Go over all segments and if it is similar, add it to the corresponding cluster
            for j = i:size(segIndex_curr,1)
                [clusters,alreadyClustered] = FnDecideSimilar(i,j,clusters,alreadyClustered, cell2mat(combined_metric(k)), thres);
            end
        end

        % Convert representation from cluster_id -> segment_id to segment_id -> cluster_id
        cluster_curr = alreadyClustered(:,1);
        cluster_segs(k,1) = {cluster_curr};

        % Add local clusters to trace data structure 
        for i=1:length(trace)
            chpoints_var = cell2mat(trace(i).chpoints_per_var(k,1));
            len_segs = length(chpoints_var)-1;
            trace(i).labels_trace_per_var(k,1) = {cluster_curr(1:len_segs,1)};
            cluster_curr(1:len_segs,:) = [];
        end
    end
end

function [cluster_global, trace] = computeClustersGlobal(x,trace, cluster_segs, segIndex, segIndex_var)
% computeClustersGlobal returns global clusters by merging local clusters
% to global clusters by considering all occuring combinations of local cluster ids.
    global num_var useLMIrefine
    % Merge combination of local cluster id to single global cluster id by
    % utilizing a map
    indices = ones(num_var,1);
    nextid = 1;
    cluster_global = zeros(size(segIndex,1),1);
    M = containers.Map('KeyType','char','ValueType','double');
    % Iterate over all present global segments
    for i=1:size(segIndex,1)
        % Create key of current global segments for map as concatenation of 
        % all corresponding local cluster ids seperated by dashs
        key = '';
        for k = 1:num_var
            if(k ~= 1)
                key = [key '-'];
            end
            cluster_curr = cell2mat(cluster_segs(k,1));
            key = [key num2str(cluster_curr(indices(k,1),1))];
            seg_curr = cell2mat(segIndex_var(k,1));
            % If for next global segment the next local segment needs to be
            % considered, update the corresponding local index
            if(segIndex(i,2) == seg_curr(indices(k,1),2))
                indices(k,1) = indices(k,1) + 1;
            end
        end
        % If key is not part of map, add the key to the map
        if(isKey(M,key) == 0)
            M(key) = nextid;
            nextid = nextid + 1;
        end
        % Save global cluster id of currently considered global cluster
        cluster_global(i,1) = M(key);
    end

    %Use LMI to possibly merge incorrectly split-up clusters
    if (useLMIrefine)
        cluster_global = refineClustersLMI(x,cluster_global,segIndex);
    end
end