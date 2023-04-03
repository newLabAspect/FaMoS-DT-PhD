function [X,Y,states] = FnTraceToTrainingData(trace, N)
%FNTRACETOTRAININGDATA creates training data for DTL out of analyzed trace
%   Now uses the fixed interval length approach to create said data
    global fixedIntervalLength precisionDTL num_var useTime
    %% sample sequence of states (fixed interval approach)
    states = [];
    values = [];
    timeSwitch = [];
    lastswitch = 1;
    indxStates = 1;
    indxstep = fixedIntervalLength;
    indx = 1;
    while true
        %finished if end is reached
        if(indx >= size(trace.x,1))
            break;
        end
        %if indx is in new segment, update indxStates (of trace.chpoints)
        if(indxStates+1 <= size(trace.chpoints,1) && indx >= trace.chpoints(indxStates+1,1))
            indxStates = indxStates + 1;
            lastswitch = indx;
        end
        %append state of current indx to states
        states = [states; trace.labels_trace(indxStates,1)];
        %interestingly if end of fixed intervals used, worse results
        %possibly because that indx might already be in next state?
        values = [values; FnRoundToInterval(trace.x(indx,1:num_var),precisionDTL)];
        timeSwitch = [timeSwitch; indx-lastswitch];
        %increase indx (of trace.x)
        indx = indx + indxstep;
    end
    %% create training data for DTL
    X = [];
    Y = [];
    %i represent all possible start positions
    for i = 1:(size(states,1)-N)
        %input to DTL are the N next states starting at i
        %append also last datapoint to X (better learning transitions)
        if (useTime)
            X = [X; states(i:(i+N-1),1)', values(i+N-1,:), timeSwitch(i+N-1,:)];
        else
            X = [X; states(i:(i+N-1),1)', values(i+N-1,:)];
        end
        %wanted output to DTL is the N'th state after i
        Y = [Y; states(i+N,1)'];
    end 
end

