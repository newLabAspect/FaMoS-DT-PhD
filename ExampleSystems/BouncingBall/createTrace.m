
for curr = 1:5

    sim("sf_bounce.slx");
    
    xout = simout;
    chpoints = [1, length(simout)];
    states = [1];
    u
    save(['ExampleSystems\BouncingBall\training', int2str(curr),'.mat'],'xout','states','chpoints');
end