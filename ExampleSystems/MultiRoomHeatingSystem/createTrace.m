Ts = 0.01;
x_init = [20.0,24.0,23.0; 21.0,23.0,25.0; 19.25,22.5,21.0; 23.5,20.0,25.75; 23.0,24.5,21.0;
          18.0,22.0,24.0; 21.0,18.5,25.0; 22.0,24.0,19.0; 21.0,25.0,22.0; 20.0,24.0,22.0];
states_init = [0;0;3;2;3;1;2;3;0;0];

h_on = 0.2; c_on = 0.3; h_off = -0.1; c_off = 0.0;
t_low_1 = 19; t_low_2 = 20; t_low_3 = 21;
t_high_1 = 24; t_high_2 = 25; t_high_3 = 26;
g_1 = 17; g_2 = 18; g_3 = 19;
d_1 = 3; d_2 = 3; d_3 = 3;

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
states = [state];
chpoints = [i];
x = [x_init(curr, :), 0.0, 0.0, 0.0];

while t < 20.0
    %dynamic
    if state == 0
        x_dot_1 = h_off * x(i,1) + c_off;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = h_off * x(i,2) + c_off;
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = h_off * x(i,3) + c_off;
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 1
        x_dot_1 = h_on * x(i,1) + c_on;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = h_off * x(i,2) + c_off;
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = h_off * x(i,3) + c_off;
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 2
        x_dot_1 = h_off * x(i,1) + c_off;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = h_on * x(i,2) + c_on;
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = h_off * x(i,3) + c_off;
        x_3 = x(i,3) + x_dot_3*Ts;
    elseif state == 3
        x_dot_1 = h_off * x(i,1) + c_off;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = h_off * x(i,2) + c_off;
        x_2 = x(i,2) + x_dot_2*Ts;
        x_dot_3 = h_on * x(i,3) + c_on;
        x_3 = x(i,3) + x_dot_3*Ts;
    end

    %transitions
    state_old = state;
    if state == 0
        if x_3 <= t_low_3
            state = 3;
        elseif x_2 <= t_low_2
            state = 2;
        elseif x_1 <= t_low_1
            state = 1;
        end
    elseif state == 1
        if x_1 >= t_high_1
            state = 0;
        elseif x_2 <= g_2 && x_1 - x_2 < d_2 
            state = 2;
        elseif x_3 <= g_3 && x_1 - x_3 < d_3 
            state = 3;
        end
    elseif state == 2
        if x_2 >= t_high_2
            state = 0;
        elseif x_1 <= g_1 && x_2 - x_1 < d_1 
            state = 1;
        elseif x_3 <= g_3 && x_2 - x_3 < d_3 
            state = 3;
        end
    elseif state == 3
        if x_3 >= t_high_3
            state = 0;
        elseif x_1 <= g_1 && x_3 - x_1 < d_1 
            state = 1;
        elseif x_2 <= g_2 && x_3 - x_2 < d_2 
            state = 2;
        end
    end

    if(state ~= state_old)
        states = [states; state];
        chpoints = [chpoints; i];
    end

    x = [x; x_1, x_2, x_3, x_dot_1, x_dot_2, x_dot_3];
    i = i +1;
    t = t + Ts;
end

xout = x(:,1:3);
chpoints = [chpoints; i];
save(['ExampleSystems\MultiRoomHeatingSystem\training', int2str(curr),'.mat'],'xout','states','chpoints');
if(curr >= 11)
    hold on;
    subplot(3,1,1);
    plot(x(:,1));
    subplot(3,1,2);
    plot(x(:,2));
    subplot(3,1,3);
    plot(x(:,3));
    hold off;
end
end