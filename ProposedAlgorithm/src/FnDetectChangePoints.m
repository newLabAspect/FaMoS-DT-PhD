function trace = FnDetectChangePoints(xout, udout, xout_shifts)
% FnDetectChangePoints returns local and global changepoints present in
% input and output traces
%   A sliding window approach, utilizing the euclidean distance between
%   the immediate past and immediate future to detect changes in dynamic
%   behavior on all available derivatives, is used
    global num_var num_ud max_deriv chp_depths
    warning('off','signal:findpeaks:largeMinPeakHeight');

    chpoints = []; % Global changepoints
    chp_depths = zeros(max_deriv+1,1); % Purely for debugging
    chp_var = cell(num_var+num_ud,1); % Local changepoints (per variable)
    for i = 1:num_var
        % Find changepoints on i'th trace incl. derivatives
        new_chp = findChangePoints(xout(:,i:num_var:end),0,1,size(xout,1),max_deriv);
        % Add newly discovered changepoints to global set of changepoints
        chpoints = union(chpoints, new_chp);
        % Add newly discovered changepoints to local set of changepoints
        % (one set per output variable)
        chp_var(i) = {sort(new_chp)};
    end
    % Repeat for input variables but only on original trace
    for i = 1:num_ud
        new_chp = findChangePoints(udout(:,i),0,1,size(udout,1),0);
        chpoints = union(chpoints, new_chp);
        chp_var(num_var+i) = {new_chp};
    end
    
    % Filter changepoints given close proximity while keeping entries
    % in global and local changepoint set consistent
    [xout,udout,xout_shifts,chpoints,chp_var] = filterChangePoints(xout, udout, xout_shifts, chpoints, chp_var);

    % Create trace datastructure
    trace.x = xout;
    trace.xs = xout_shifts;
    trace.chpoints = chpoints;
    trace.chpoints_per_var = chp_var;
    trace.ud = udout;
    trace.labels_num = []; 
    trace.labels_trace = [];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist = computeDistance(der)
% computeDistance returns the similarity of the immediate past to the
% immediate future for all inner indices using the euclidean distance metric
    global windowSize
    dist = zeros(windowSize,1);
    % Iterate over all inner indices usable in sliding window approach
    for i = (windowSize+1):(length(der)-windowSize)
        before = der((i-windowSize):(i-1));
        after = der((i+1):(i+windowSize));
        % Calculate euclidean distance between immediate past and future (offsets are removed)
        dist_new = sum(abs(before-before(1)-(after-after(1))));
        dist = [dist; dist_new];
    end
end

function locs = findChangePoints(xout,depth,starting,ending,max_depth)
% findChangePoints returns changepoints present in an interval using all
% available derivatives up to the selected one
%   The search is realized recursively: Terminal cases are if the interval
%   is too small or if the selected derivative is not avialable. Otherwise,
%   changepoints are detected on the selected derivative and used to
%   formulate new intervals which are examined in the next recursion with
%   the next derivative. All collected changepoints are merged to one set.
    global windowSize chp_depths

    % Check for terminal cases 
    locs = [];
    if depth > max_depth || ending-starting-1 < 2*windowSize
        return;
    end

    der = xout(starting:ending,depth+1);
    dist = computeDistance(der);

    % Peaks in distance indicate change in dynamic behavior (introduce para for minimal peak height?)
    [~,locsDist] = findpeaks(dist,"MinPeakHeight",5);
    % Convert index on interval to index on trace
    locsHere = locsDist+starting-1;
    locsHere = sort(locsHere);
    % Remove changepoints that likely come from same system mode switch (why 1.5*... ?)
    locsHere = filterindx(locsHere,1.5*windowSize);
    chp_depths(depth+1) = chp_depths(depth+1) + length(locsHere); % Debugging
    locs = [locs; locsHere];

    % Formulate new intervals for next recursive calls and add changepoints
    % found by them to set changepoints
    locsHere = [starting-windowSize/2; locsHere; ending+windowSize/2];
    for i = 1:(length(locsHere)-1)
        newStart = (locsHere(i)+windowSize/2);
        newEnd = (locsHere(i+1)-windowSize/2);
        locsNew = findChangePoints(xout,depth+1,newStart,newEnd,max_depth);
        if(size(locsNew,1) ~= 0)
            locs = [locs; locsNew(:,1)];
        end
    end

    % Start and end of trace are also changepoints, i.e., add them to set
    % of changepoints if this is the top-most function call
    if depth == 0
        locs = [1; locs; length(der)];
    end
end

function [xout,udout,xout_shifts,chpoints,chp_var] = filterChangePoints(xout, udout, xout_shifts, chpoints, chp_var)
% filterChangePoints returns xout, udout and the local and global 
% changepoint set after applying some filtering rules
%   Changepoints that are visible on mutliple output variables and thus 
%   generate multiple entries in the global changepoint set are merged;
%   short last segments are deleted; entries in local changepoint sets
%   changed so that positions are consistent with the global set
    global windowSize num_var num_ud
    % Remove changepoints in global set that were detected on multiple output variables
    chpoints = filterindx(chpoints,windowSize);

    % Remove last segment if too short
    if(chpoints(end,1)-chpoints(end-1,1) < 2 * windowSize)
        xout = xout(1:chpoints(end-1,1),:);
        xout_shifts = xout_shifts(1:chpoints(end-1,1),:);
        if num_ud ~= 0
            udout = udout(1:chpoints(end-1,1),:);
        end
        chpoints = chpoints(1:(end-1),1);
        for i = 1:(num_var+num_ud)
            current_chps = cell2mat(chp_var(i));
            % Remove last entry (the same for all output variables)
            current_chps = current_chps(1:(end-1),1);
            % Add new common end point (potential duplicate removed in next step)
            current_chps = [current_chps; chpoints(end,1)];
            chp_var(i) = {current_chps};
        end
    end

    % Ensure that if multiple changepoints in global set are merged, the
    % same unique changepoint value will be used in local sets
    for i = 1:(num_var+num_ud)
        current_chps = cell2mat(chp_var(i));
        % Find closest changepoint in global set (distance can be zero) and
        % use this changepoint value going forward
        for j = 1:length(current_chps)
            [~,ind] = min(abs(chpoints-current_chps(j,1)));
            current_chps(j,1) = chpoints(ind,1);
        end
        % Remove possible duplicates from previous step
        current_chps = unique(current_chps);
        chp_var(i) = {current_chps};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Carried over from HAutLearn (with one minor bugfix)
function indx = filterindx(indx,windw)
% filterindx checks if there are two changepoint inside a windw and erases 
% the 2nd one if this is the case
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

