Ts = 0.01;
x_init = [15.0,10.0; 10.0,10.0; 23.0,14.0; 18.0,13.0; 16.0,12.0;
          24.0,13.0; 19.0,12.0; 12.0,10.0; 22.5,13.5; 19.5,11.0];
states_init = [1; 1; 2; 3; 1; 2; 3; 1; 2; 3];

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
x = [x_init(curr, :), 0.0, 0.0];

while t < 20.0
    if state == 1
        x_dot_1 = 0.5 * x(i,1) + 0.5;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = 0.3 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
    elseif state == 2
        x_dot_1 = - 0.5 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = 0.3 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
    elseif state == 3
        x_dot_1 = - 0.5 * x(i,1);
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = - 0.5 * x(i,2);
        x_2 = x(i,2) + x_dot_2*Ts;
    end

    if x_1 >= 25 && state == 1
        state = 2;
    elseif x_1 <= 20 && state == 2
        state = 3;
    elseif x_2 <= 10 && state == 3
        state = 1;
    end

    x = [x; x_1, x_2, x_dot_1, x_dot_2];
    i = i +1;
    t = t + Ts;
end

xout = x(:,1:2);
save(['evaluation\MultiVariablesHeatingSystem\training', int2str(curr),'.mat'],'xout');

end