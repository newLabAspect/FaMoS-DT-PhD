Ts = 0.01;
x_init = [10.0,8.0,0.0; 14.0,12.0,0.0; 18.0,-6.0,0.0; 20.0,6.0,0.0; 22.0,2.0,0.0;
          23.5,8.0,0.0; 22.0,-4.0,0.0; 24.0,6.0,0.0; 26.0,-8.0,0.0; 23.0,-2.0,0.0];
states_init = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
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
    end

    if x_ >= 25 && state ==1
        state = 2;
    elseif x_ <= 15 && state == 2
        state = 1;
    end

    x = [x; x_, x_dot, x_dot_dot];
    i = i +1;
    t = t + Ts;
end

xout = x(:,1);
save(['evaluation\TwoStateHybridAutomaton\training', int2str(curr),'.mat'],'xout');

end