Ts = 0.01;
x_init = [10.0,8.0,0.0; 14.0,12.0,0.0; 18.0,-6.0,0.0; 20.0,6.0,0.0; 20.0,0.0,0.0;
          23.5,8.0,0.0; 23.0,-2.0,0.0; 24.0,6.0,0.0; 19.0,0.0,0.0; 18.0,0.0,0.0];
states_init = [1; 1; 1; 1; 3; 2; 2; 2; 3; 3];

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
states = [state];
chpoints = [i];
x = x_init(curr, :);

while t < 20.0
    if state == 1
        x_dot_dot = 0.5*x(i,1) - 0.5*x(i,2);
        x_dot = x(i,2) + x_dot_dot*Ts;
        x_ = x(i,1) + x_dot*Ts;
    elseif state == 2
        x_dot_dot = -0.5*x(i,1) - 0.5*x(i,2);
        x_dot = x(i,2) + x_dot_dot*Ts;
        x_ = x(i,1) + x_dot*Ts;
    elseif state == 3
        x_dot = -0.5 * x(i,1);
        x_ = x(i,1) + x_dot*Ts;
    end

    state_old = state;
    if x_ >= 25 && state ==1
        state = 2;
    elseif x_ <= 20 && state == 2
        state = 3;
    elseif x_ <= 15 && state == 3
        state = 1;
        x_dot_dot = (x(i,2)-x(i-5,2))/(5*Ts);
    end

    if(state ~= state_old)
        states = [states; state];
        chpoints = [chpoints; i];
    end

    x = [x; x_, x_dot, x_dot_dot];
    i = i +1;
    t = t + Ts;
end

xout = x(:,1);
chpoints = [chpoints; i];
save(['ExampleSystems\ThreeStateHybridAutomaton\training', int2str(curr),'.mat'],'xout','states','chpoints');

end