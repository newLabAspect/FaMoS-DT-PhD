Ts = 0.01;
x_init = [23,25.0,28.0; 22.0,21.0,25.0; 23.5,20.0,14.0; 18.0,13.0,21.0; 18.0,8.0,13.0;
          14.0,18.0,17.0; 6.0,22.0,14.0; 12.0,8.0,23.0; 4.0,2.0,4.0; 0.5,3.0,1.0];
states_init = [0; 0; 1; 2; 3; 4; 5; 6; 7; 7];

l_1 = -0.08; l_2 = -0.07; l_3 = -0.05;
f_1 = 11; f_2 = 10; f_3 = 9;
t_f = 20; t_h = 15; t_m = 10; t_l = 5;
states = cell(10,1);

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
states(curr) = {state};
allstates = [state];
x = [x_init(curr, :), 0.0, 0.0, 0.0];

while t < 20.0
    %dynamic
    if state == 0
        x_dot_1 = l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 1
        x_dot_1 = l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = f_3 + l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 2
        x_dot_1 = l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = f_2 + l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 3
        x_dot_1 = l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = f_2 + l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = f_3 + l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 4
        x_dot_1 = f_1 + l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 5
        x_dot_1 = f_1 + l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = f_3 + l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 6
        x_dot_1 = f_1 + l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = f_2 + l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 7
        x_dot_1 = f_1 + l_1 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = f_2 + l_2 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = f_3 + l_3 * x(i,3);
        x_3 = x(i,3) + x_dot_3*Ts;
    end

    %transitions
    state_old = state;
    if state == 0
        if x_3 <= t_h
            state = 1;
        elseif x_2 <= t_h
            state = 2;
        elseif x_1 <= t_h
            state = 4;
        end
    elseif state == 1
        if x_3 >= t_f
            state = 0;
        elseif x_2 <= t_m
            state = 3;
        elseif x_1 <= t_m
            state = 5;
        end
    elseif state == 2
        if x_3 <= t_m
            state = 3;
        elseif x_2 >= t_f
            state = 0;
        elseif x_1 <= t_m
            state = 6;
        end
    elseif state == 3
        if x_3 >= t_f
            state = 2;
        elseif x_2 >= t_f
            state = 1;
        elseif x_1 <= t_l
            state = 7;
        end
    elseif state == 4
        if x_3 <= t_m
            state = 5;
        elseif x_2 <= t_m
            state = 6;
        elseif x_1 >= t_f
            state = 0;
        end
    elseif state == 5
        if x_3 >= t_f
            state = 4;
        elseif x_2 <= t_l
            state = 7;
        elseif x_1 >= t_f
            state = 1;
        end
    elseif state == 6
        if x_3 <= t_l
            state = 7;
        elseif x_2 >= t_f
            state = 4;
        elseif x_1 >= t_f
            state = 2;
        end
    elseif state == 7
        if x_3 >= t_f
            state = 6;
        elseif x_2 >= t_f
            state = 5;
        elseif x_1 >= t_f
            state = 3;
        end
    end

    if(state ~= state_old)
        states(curr) = {[cell2mat(states(curr)), state]};
    end
    allstates = [allstates, state];

    x = [x; x_1, x_2, x_3, x_dot_1, x_dot_2, x_dot_3];
    i = i +1;
    t = t + Ts;
end

xout = x(:,1:3);
save(['evaluation\ComplexTankSystem\training', int2str(curr),'.mat'],'xout');
hold on;
subplot(4,1,1);
plot(allstates);
subplot(4,1,2);
plot(x(:,1));
subplot(4,1,3);
plot(x(:,2));
subplot(4,1,4);
plot(x(:,3));
hold off;
end