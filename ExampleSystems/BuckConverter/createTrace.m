Ts = 1E-5;
%x_1 is il while x_2 is v_c
x_init = [2.0,7.0; 8.0,2.0; 14.0,8.0; 20.0,14.0; 26.0,12.0;
          -0.05,12.5; -0.05,14.0; -0.05,16; 1.0,8.0; 4.0,4.0];
states_init = [1;1;1;2;2;3;3;3;1;1];

a00c = -271.6981; a01c = -377.3585; a10c = 454.5455; a11c = -45.4545; b0c = 377.3585; b1c = 0;
a00o = -196.2264; a01o = -377.3585; a10o = 454.5455; a11o = -45.4545; b0o = 0; b1o = 0;
Vs = 24; VcH = 12.1; VcL = 11.9;

for curr = 1:length(states_init)

i = 1;
t = 0.0;
state = states_init(curr,1);
states = [state];
chpoints = [i];
x = [x_init(curr, :), 0.0, 0.0];

while t < 0.02
    if state == 1
        x_dot_1 = a00c * x(i,1) + a01c * x(i,2) + b0c * Vs;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = a10c * x(i,1) + a11c * x(i,2) + b1c * Vs;
        x_2 = x(i,2) + x_dot_2*Ts;
    elseif state == 2
        x_dot_1 = a00o * x(i,1) + a01o * x(i,2) + b0o * Vs;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = a10o * x(i,1) + a11o * x(i,2) + b1o * Vs;
        x_2 = x(i,2) + x_dot_2*Ts;
    elseif state == 3
        x_dot_1 = 0.0;
        x_1 = x(i,1) + x_dot_1*Ts;
        x_dot_2 = a11o * x(i,2) + b1o * Vs;
        x_2 = x(i,2) + x_dot_2*Ts;
    end

    state_old = state;
    if x_2 >= VcH && state == 1
        state = 2;
    elseif x_1 <= 0 && state == 2
        state = 3;
    elseif x_2 <= VcL && (state == 2 || state == 3)
        state = 1;
    end

    if(state ~= state_old)
        states = [states; state];
        chpoints = [chpoints; i];
    end

    x = [x; x_1, x_2, x_dot_1, x_dot_2];
    i = i +1;
    t = t + Ts;
end

xout = x(:,1:2);
chpoints = [chpoints; i];
save(['ExampleSystems\BuckConverter\training', int2str(curr),'.mat'],'xout','states','chpoints');

end