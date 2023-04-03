function trace = FnDetectChangePoints(xout, num_var)
%FNDETECTCHANGEPOINTS Returns change points of xout
%   Now uses the more sophisticated DTW approach to find changepoints
    global windowSize max_deriv chp_depths
    warning('off','signal:findpeaks:largeMinPeakHeight');

    chpoints = [];
    chp_depths = zeros(max_deriv+1,1);
    for i = 1:num_var
        %add new found changepoints of i'th trace to set
        chpoints = union(chpoints, findChangePoints(xout(:,i:num_var:end),0,1,size(xout,1)));
    end
    
    % remove redundant chpoints
    chpoints = filterindx(chpoints,windowSize);
    %xout_reduced= xout(:, 1:num_var);

    %remove last segment if too short
    if(chpoints(end,1)-chpoints(end-1,1) < 2 * windowSize)
        xout = xout(1:chpoints(end-1,1),:);
        chpoints = chpoints(1:end-1,1);
    end

    trace.x = xout;
    trace.chpoints = chpoints;
    trace.ud = [];
    trace.labels_num = []; 
    trace.labels_trace = [];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function locs = findChangePoints(xout,depth,starting,ending)
    global Ts windowSize max_deriv chp_depths
    locs = [];
    if depth > max_deriv || ending-starting-1 < 2*windowSize
        return;
    end
    der = xout(starting:ending,depth+1);

    %distanceMetric of DTW
    dist = zeros(windowSize,1);
    %diagonility metric of DTW
    diag = zeros(windowSize,1);
    %because before and after points are needed keep distance to bounds
    for i = (windowSize+1):(length(der)-windowSize)
        %points right before i
        before = der((i-windowSize):(i-1));
        %points right after i
        after = der((i+1):(i+windowSize));
        %calc DTW (corrected by initial values)
        [dist_new, i_x, i_y] = dtw(before-before(1),after-after(1));
        %update distance metric
        dist = [dist; dist_new];
        %update diagnoality metric
        diag_new = corrcoef(i_x,i_y);
        diag = [diag; diag_new(1,2)];
    end

    %find dips in diagonality (maybe move minpaekheight to paras?)
    %[pksDiag,locsDiag] = findpeaks(1-diag,"MinPeakHeight",0.05);
    %find peaks in distance (maybe move mindipdepth to paras?)
    [pksDist,locsDist] = findpeaks(dist,"MinPeakHeight",5);
    %both are possible changepoints
    locsHere = locsDist+starting-1; %union(locsDiag,locsDist);
    locsHere = sort(locsHere);
    %get rid of too close changepoints
    locsHere = filterindx(locsHere,1.5*windowSize);
    chp_depths(depth+1) = chp_depths(depth+1) + length(locsHere);
    locs = [locs; locsHere]; %, depth*ones(size(locsHere))
    locsHere = [starting-windowSize/2; locsHere; ending+windowSize/2];
    for i = 1:(length(locsHere)-1)
        newStart = (locsHere(i)+windowSize/2);
        newEnd = (locsHere(i+1)-windowSize/2);
        locsNew = findChangePoints(xout,depth+1,newStart,newEnd);
        if(size(locsNew,1) ~= 0)
            locs = [locs; locsNew(:,1)]; %, locsNew(:,2)
        end
    end

    if depth == 0
        locs = union(1,[locs; length(der)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indx = filterindx(indx,windw)
%if there are two changepoint inside a windw the 2nd one gets erased
    n = 1;
    while true
        if n >= length(indx)
            break;
        end
        id1 = indx(n);
        while true
            if n+1 > length(indx) %not >= bc otherwise last point not checked
                break;
            end
            id2 = indx(n+1);
            if id2-id1<=windw
                indx(n+1) = [];
            else
                break;
            end
        end
        n = n+1;
    end
end

