function [X,Y,states] = FnTraceToTrainingData(trace)
%FnTraceToTrainingData generates feature vectors and class labes from a
%clustered trace object
%   The feature vector includes all trace values, the state at the current 
%   position and possibly the time since the last system mode switch while 
%   the label includes the state at the next position
    global precisionDTL num_var useTime

    % Preallocate matrices needed to save samples
    states = zeros(size(trace.x,1)-1,1);
    values = zeros(size(trace.x,1)-1,num_var);
    timeSwitch = zeros(size(trace.x,1)-1,1);
    % Variables needed to keep track of states
    lastswitch = 1;
    indxStates = 1;

    for indx = 1:(size(trace.x,1)-1)
        % Current index associated with system mode switch, update corresponding vars
        if(indxStates+1 <= size(trace.chpoints,1) && indx >= trace.chpoints(indxStates+1,1))
            indxStates = indxStates + 1;
            lastswitch = indx;
        end

        % Save samples needed for feature vector and class label creation
        states(indx,:) = [trace.labels_trace(indxStates,1)];
        values(indx,:) = [trace.x(indx,1:num_var)];
        timeSwitch(indx,:) = indx-lastswitch;
    end

    % Create matrices containing feature vectors and corresponding class labels
    points = 1:(size(states,1)-1);
    if (useTime)
        X = [states(points,1), values(points,:), timeSwitch(points,1)];
    else
        X = [states(points,1), values(points,:)];
    end
    Y = [states(points+1,1)];
end