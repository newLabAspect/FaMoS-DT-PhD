correct = 0;
false = 0;

for i = [2,6]
    for j = 2:(length(trace(i).chpoints)-1)
        x1 = trace(i).x(trace(i).chpoints(j),1);
        x2 = trace(i).x(trace(i).chpoints(j),2);
        x3 = trace(i).x(trace(i).chpoints(j),3);
        state = trace(i).labels_trace(j-1);
        state_pred = state;
        state_next = trace(i).labels_trace(j);

        %conditions are extracted from automata_learning.xml

        % 8 trace used
        if(state == 2 && -1.159 * x1 - 0.007756 * x2 - 0.007662 * x3 + 1.0 < 0.01)
            state_pred = 1;
        elseif(state == 1 && 0.0008079 * x3 - 1.667 * x2 - 0.0005141 * x1 + 1.0 > -0.01)
            state_pred = 4;
        elseif(state == 4 && 0.002552 * x1 - 1.248 * x2 - 0.001274 * x3 + 1.0 < 0.01)
            state_pred = 1;
        elseif(state == 1 && 0.001704 * x2 - 1.568 * x1 - 0.0008483 * x3 + 1.0 > -0.01)
            state_pred = 2;
        elseif(state == 3 && 0.0002276 * x1 + 0.0006757 * x2 - 1.399 * x3 + 1.0 < 0.01)
            state_pred = 1;
        elseif(state == 1 && -0.000688 * x1 - 0.0002388 * x2 - 1.866 * x3 + 1.0 > -0.01)
            state_pred = 3;
        elseif(state == 4 && 0.002239 * x3 - 1.216 * x2 - 0.0446 * x1 + 1.0 < 0.01)
            state_pred = 2;
        elseif(state == 2 && -0.8567 * x1 - 0.04278 * x2 - 0.4475 * x3 + 1.0 < 0.01)
            state_pred = 3;
        end

        if(state_pred == state_next)
            correct = correct + 1;
        else
            false = false +1;
        end
    end
end