Ts = 0.01;
x_init = [10.0,8.0,0.0; 14.0,12.0,0.0; 18.0,-6.0,0.0; 20.0,6.0,0.0; 20.0,0.0,0.0;
          23.5,8.0,0.0; 23.0,-2.0,0.0; 24.0,6.0,0.0; 19.0,0.0,0.0; 18.0,0.0,0.0];
states_init = [1; 1; 1; 1; 3; 2; 2; 2; 3; 3];

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
states = [state];
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

    if x_ >= 25 && state == 1
        toss = binornd(1,0.5);
        if toss == 1 % success
            state = 3;
        else
            state = 2;
        end
        states = [states, state];
    elseif x_ <= 20 && state == 2
        state = 1;
        states = [states, state];
    elseif x_ <= 15 && state == 3
        state = 1;
        states = [states, state];
        x_dot_dot = (x(i,2)-x(i-5,2))/(5*Ts);
    end

    x = [x; x_, x_dot, x_dot_dot];
    i = i + 1;
    t = t + Ts;
end

xout = x(:,1);
disp(states);
%to prevent accidents
%save(['evaluation\ThreeStateStochasticHybridAutomaton\training', int2str(curr),'.mat'],'xout');

end