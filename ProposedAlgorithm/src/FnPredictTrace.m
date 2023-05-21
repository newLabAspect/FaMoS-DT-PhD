function [sim_x, sim_state] = FnPredictTrace(trace,Mdl,ode)
    global fixedIntervalLength useTime num_var Ts precisionDTL
    offsetPred = 0+1;
    sim_x = trace.x(offsetPred,:);
    lastSwitch = 1;
    sim_state = 1;
    curr_state = 1;
    for j = 1:size(trace.chpoints,1)
        if trace.chpoints(j,1) > offsetPred
            break;
        end
        lastSwitch = trace.chpoints(j,1);
        sim_state = trace.labels_trace(j,1);
        curr_state = trace.labels_trace(j,1);
    end
    for i = ((offsetPred-1)/fixedIntervalLength+1):floor(size(trace.x)/fixedIntervalLength)
        A = cell2mat(ode(curr_state));
        A = A(1:size(A,1),1:size(A,1));
        curr_x = sim_x(end,1:size(A,1))';
        new_x_dot = A * curr_x;
        new_x = [curr_x(1:end,1); new_x_dot((end-num_var+1):end,1)];
        for j = (size(new_x,1)-num_var):-1:1
            new_x(j,1) = new_x(j,1) + new_x(j+num_var,1) * Ts * fixedIntervalLength;
        end
        sim_x(end,size(A,1)+(1:num_var)) = new_x_dot((end-num_var+1):end,1)';
        sim_x = [sim_x; new_x', zeros(1,size(sim_x,2)-size(new_x,1))];
        sim_state = [sim_state; curr_state];
        last_state = curr_state;
        if (useTime)
            curr_state = predict(Mdl,[last_state FnRoundToInterval(new_x(1:num_var,1)',precisionDTL) (i-1)*fixedIntervalLength-lastSwitch]); 
        else
            curr_state = predict(Mdl,[last_state FnRoundToInterval(new_x(1:num_var,1)',precisionDTL)]);
        end
        if(last_state ~= curr_state)
            lastSwitch = i;
        end
    end
end