function trace = FnClusterSegs(trace, x, ud)
    global winlen thresClusterMax thresClusterMin num_var offsetCluster facThres;

    segIndex = [0,0];
    segIndex_var = cell(num_var,1);
    segIndex_var(1:num_var,1) = {[0,0]};
    for i=1:length(trace)  
        %segments globally
        chpoints = (trace(i).chpoints);
        %start and end point of each segment are saved
        chsegments = [chpoints(1:end-1), chpoints(2:end)];
        chsegments(2:end,1) = chsegments(2:end,1)+1;
        %in x all traces are appended thus offsets need to be added for
        %changepoints to be consistent with x (and ud)
        segIndex = [segIndex; segIndex(end,2) + chsegments];

        %segments per output variable
        chp_var = (trace(i).chpoints_per_var);
        for j=1:num_var
            %start and end position of each segment saved
            chp_curr = cell2mat(chp_var(j,1));
            chsegments_var = [chp_curr(1:end-1),chp_curr(2:end)];
            chsegments_var(2:end,1) = chsegments_var(2:end,1)+1;
            %offset so that indices are relative to x
            segIndex_var_temp = cell2mat(segIndex_var(j,1));
            segIndex_var(j,1) = {[segIndex_var_temp; segIndex_var_temp(end,2) + chsegments_var]};
        end 
    end
    %dummy line is deleted (used so that first append offset is 0)
    segIndex(1,:) = [];
    for j=1:num_var
        segIndex_var_temp = cell2mat(segIndex_var(j,1));
        segIndex_var_temp(1,:) = [];
        segIndex_var(j,1) = {segIndex_var_temp};
    end
    
    %create metric
    dist_all = cell(num_var,1);
    diag_all = cell(num_var,1);
    combined_metric = cell(num_var,1);
    for k = 1:num_var
        segIndex_curr = cell2mat(segIndex_var(k,1));
        dist_curr = zeros(size(segIndex_curr,1),size(segIndex_curr,1));
        diag_curr = zeros(size(segIndex_curr,1),size(segIndex_curr,1));
        combined_curr = zeros(size(segIndex_curr,1),size(segIndex_curr,1));
        for i = 1:size(segIndex_curr,1)
            start_i = segIndex_curr(i,1)+winlen;
            end_i = segIndex_curr(i,2)-winlen;
            seg_i = x(start_i:end_i,k+offsetCluster);
            for j = i:size(segIndex_curr,1)
                start_j = segIndex_curr(j,1)+winlen;
                end_j = segIndex_curr(j,2)-winlen;
                seg_j = x(start_j:end_j,k+offsetCluster);
                common_len = min(size(seg_i,1),size(seg_j,1));
                
                incConst = 1;
                if offsetCluster ~= 0
                    incConst = 0;
                end
                
                [dist_end, i_x, i_y] = dtw(seg_i((end-common_len+1):end,1)-seg_i(end,1)*incConst,seg_j((end-common_len+1):end,1)-seg_j(end,1)*incConst);
                diag = corrcoef(i_x,i_y);
                diag_end = diag(1,2);
                dist_end = dist_end/size(i_x,2);
                %diag can be NaN when sigma_x or sigma_y is 0, if this is the
                %case diag can be set to 0 bc one const. is repeated in i_.
                if isnan(diag_end)
                    diag_end = 0.0;
                end
                sim_index_end = 0.5*dist_end+0.5*(1-diag_end);

                [dist_start, i_x, i_y] = dtw(seg_i(1:common_len,1)-seg_i(1,1)*incConst,seg_j(1:common_len,1)-seg_j(1,1)*incConst);
                diag = corrcoef(i_x,i_y);
                diag_start = diag(1,2);
                dist_start = dist_start/size(i_x,2);
                %diag can be NaN when sigma_x or sigma_y is 0, if this is the
                %case diag can be set to 0 bc one const. is repeated in i_.
                if isnan(diag_start)
                    diag_start = 0.0;
                end
                sim_index_start = 0.5*dist_start+0.5*(1-diag_start);

                if (sim_index_start <= sim_index_end)
                    dist_curr(i,j) = dist_start;
                    diag_curr(i,j) = diag_start;
                    combined_curr(i,j) = sim_index_start;
                    %Due to symmetry
                    dist_curr(j,i) = dist_start;
                    diag_curr(j,i) = diag_start;
                    combined_curr(j,i) = sim_index_start;
                else
                    dist_curr(i,j) = dist_end;
                    diag_curr(i,j) = diag_end;
                    combined_curr(i,j) = sim_index_end;
                    %Due to symmetry
                    dist_curr(j,i) = dist_end;
                    diag_curr(j,i) = diag_end;
                    combined_curr(j,i) = sim_index_end;
                end
            end
        end
        dist_all(k,1) = {dist_curr};
        diag_all(k,1) = {diag_curr};
        combined_metric(k,1) = {combined_curr};
    end

    %actual clustering
    for i=1:length(trace)
        trace(i).labels_trace_per_var = cell(num_var,1);
    end
    cluster_segs = cell(num_var,1);
    for k = 1:num_var
        segIndex_curr = cell2mat(segIndex_var(k,1));
        clusters = cell(size(segIndex_curr,1),1);
        alreadyClustered = zeros(size(segIndex_curr,1),2);
        for i = 1:size(segIndex_curr,1)
            thres = facThres*min(combined_metric{k,1}(i,union(1:(i-1),(i+1):end)));
            if(thres < thresClusterMin)
                thres = thresClusterMin;
            elseif(thres > thresClusterMax)
                thres = thresClusterMax;
            end
            for j = i:size(segIndex_curr,1) %Possibly 1 instaed of 1
                [clusters,alreadyClustered] = FnDecideSimilar(i,j,clusters,alreadyClustered, cell2mat(combined_metric(k)), thres);
            end
        end

        %assign local cluster IDs for each output variable

        cluster_curr = zeros(size(segIndex_curr,1),1);
        for i = 1:size(segIndex_curr,1)
            cluster = cell2mat(clusters(i));
            for j =1:size(cluster,2)
                cluster_curr(cluster(j)) = i;
            end
        end
        cluster_segs(k,1) = {cluster_curr};

        for i=1:length(trace)
            chpoints_var = cell2mat(trace(i).chpoints_per_var(k,1));
            len_segs = length(chpoints_var)-1;
            %add computed labels into trace datastructure
            trace(i).labels_trace_per_var(k,1) = {cluster_curr(1:len_segs,1)};
            cluster_curr(1:len_segs,:) = [];
        end
    end

    %merge local clusters to global clusters

    indices = ones(num_var,1);
    nextid = 1;
    cluster_global = zeros(size(segIndex,1),1);
    M = containers.Map('KeyType','char','ValueType','double');
    for i=1:size(segIndex,1)
        key = '';
        for k = 1:num_var
            if(k ~= 1)
                key = [key '-'];
            end
            cluster_curr = cell2mat(cluster_segs(k,1));
            key = [key num2str(cluster_curr(indices(k,1),1))];
            seg_curr = cell2mat(segIndex_var(k,1));
            if(segIndex(i,2) == seg_curr(indices(k,1),2))
                indices(k,1) = indices(k,1) + 1;
            end
        end
        if(isKey(M,key) == 0)
            M(key) = nextid;
            nextid = nextid + 1;
        end
        cluster_global(i,1) = M(key);
    end

    labels_num = unique(cluster_global(:,1));

    for i=1:length(trace)
        chpoints = (trace(i).chpoints);
        len_segs = length(chpoints)-1;
        trace(i).labels_num = labels_num;
        %add computed labels into trace datastructure
        trace(i).labels_trace = cluster_global(1:len_segs,1);
        cluster_global(1:len_segs,:) = [];
    end
end

function [clusters,alreadyClustered] = FnDecideSimilar(i,j,clusters,alreadyClustered, combined_metric, thres)
    similar = true;
    sum_metric = 0;
    if(combined_metric(i,j) > thres)
        similar = false;
    end
    sum_metric = sum_metric + combined_metric(i,j); %relict from global approach
    if(similar && (alreadyClustered(j,1) == 0 || alreadyClustered(j,2) > sum_metric))
        if(alreadyClustered(j,1) ~= 0)
            clusters(alreadyClustered(j,1)) = {setdiff(cell2mat(clusters(alreadyClustered(j,1))),j)};
        end
        if(alreadyClustered(i,1) ~= 0)
            clusters(alreadyClustered(i,1)) = {union(cell2mat(clusters(alreadyClustered(i,1))),j)};
            alreadyClustered(j,1) = alreadyClustered(i,1);
        else
            clusters(i) = {union(cell2mat(clusters(i)),j)};
            alreadyClustered(j,1) = i;
        end
        if(i ~= j)
            alreadyClustered(j,2) = sum_metric;
        end
    end
end

function labels = FnIndxarrange(cluster_idx, index)
%Problem: cluster_idx includes all segment id's that are similar to the
%first segment of the trace it was called with. Especially the segments are
%provided sorted by length. Here the connection back is made.
    %Extract the blocks with similar segments (divided by -1)
    newindx = [];
    NOne = find(cluster_idx == -1);
    %for all entries in first block add a 1 to labels
    labels = ones(NOne(1)-1,1);
    nNOne = length(NOne);
    indxx = 1:length(index);
    
    %first block
    starti = 1;
    endi = NOne(1) - 1;
    %copy all segment ids from first block into idx 
    %Note: indxx currently identity function (but transposes)
    indx_temp = indxx(cluster_idx(starti : endi));
    idx = indx_temp;
    %remove all extracted segment ids from indxx
    indxx = setdiff(indxx, indx_temp);
    for n = 2:nNOne
        %add each entry in n'th block a n two labels
        labels = [labels; n*ones(NOne(n)- NOne(n-1)-1,1);];
        %n'th block
        starti = NOne(n -1) + 1;
        endi = NOne(n) - 1;
        %copy all segment ids from n'th block into idx while keeping order
        indx_temp = indxx(cluster_idx(starti : endi)); 
        idx = [idx, indx_temp];
        %remove all extracted segment ids from indxx
        indxx = setdiff(indxx, indx_temp);
        
    end
    %Calc back to id relative to original timestamps
    idx = index(idx);
    %create one table with ids and labels
    labels = [idx, labels];
end

%chpoints maybe not best name, segments would be better
function indx = FnRecursive(x,ud,chpoints)  
    global sigma winlen;
    
    %stop recursion if len of x reaches 0
    if size(x,1) ~= 0
        
        %setup
        xc = [];
        udc = [];
        chpc = 1;
        indx = 1;
        unindx = [];
        
        %calc number of segments and catch special case 1
        num_segments = size(chpoints,1);
        if num_segments ==1
            indx = [1;-1];
            return;
        end
        
        %variable number
        num_var = size(x,1);
        %size of ud
        num_ud = size(ud,1);
        
        %extract first segment
        index_start1 = chpoints(1,1);
        index_end1 = chpoints(1,2);
        xseg1 = x(:,index_start1: index_end1);
        %xk1: delete some points in end bc of inaccuracy of chpoints
        xk1 = xseg1(:,1:end - winlen); 
        ud1 = [];
        if num_ud ==0
            %not used: ud1 has same dimensions as xk1 but all are ones
            ud1 = ones(1,size(xk1,2));
        else
            udseg1 = ud(:,index_start1: index_end1);
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
        
        %try to find similiar segments by comparing to all other not yet matched segments    
        for j = 2:num_segments
            disp([int2str(j),'/', int2str(num_segments)])
            %extract j'th segment, xkj eq. to xk1, udj eq. ud1, xkij eq.
            %xki1, Okj eq. Ok and Okij eq. Oki
            index_startj = chpoints(j,1);
            index_endj = chpoints(j,2);
                       
            xsegj = x(:,index_startj: index_endj);
            xkj = xsegj(:,1:end - winlen);
            xkij = [];
            xkij = xsegj(:, winlen+1:end);
            
            if num_ud == 0
                udj = ones(1,size(xkj,2));
            else
                udsegj = ud(:,index_startj: index_endj);
                udj = udsegj(:,1:end - winlen);
                udj = [udj; ones(1,size(xkj,2))];
            end

            Okj = [xkj; udj];
            Okij = xkij;
            
            %append to already existing LMI
            setlmis(lmi);
            error_tj = (sigma * size(Okj,2))^2;
            for i = 1:num_var
                n = (i+(j-1)*num_var);
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
                if num_segments ==2 
                    indx = [1;2;-1];
                    return;
                end
                %append j to indx
                indx = [indx; j]; 
                lenj = length(Okj);
                len = length(Ok);
                %reset LMI to before iteration
                for i = 1:size(Okij,1)
                    n = (j-1)*num_var+i;
                    lmi = dellmi(lmi,n);
                end
            %1 and j are not similar
            else
                if num_segments ==2 
                    indx = [1;-1;1;-1];
                    return;
                end
                %reset LMI to before iteration
                for i = 1:size(Okij,1)
                    n = (j-1)*num_var+i;
                    lmi = dellmi(lmi,n);
                end
                %remember that this segment needs to be redone (xc & chpc)
                xc = [xc, xsegj]; 
                if num_ud ==0
                    udc = [];
                else
                    udc = [udc, udsegj];
                end
                if chpc == 1
                    chpc = [1,size(xsegj,2)];
                else
                    chpc = [chpc;[chpc(end)+1,size(xsegj,2)+ chpc(end)]];
                end
            end
        end
        %add -1 to indx to symbolize that all segments before are similar
        indx = [indx; -1];
        %recursive call with not yet matching segments
        indxx = FnRecursive(xc,udc,chpc);
        %add results from recursive call to general result
        indx = [indx; indxx];
         
    else
        indx = [];
    end
end