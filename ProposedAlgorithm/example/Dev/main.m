%% Intital Setup and parameter settings
clc
clear

addpath(['ProposedAlgorithm', filesep, 'example', filesep, 'Dev']);
addpath(['ProposedAlgorithm', filesep, 'src']);
global sigma num_var num_ud winlen Ts Time windowSize fixedIntervalLength max_deriv thresClusterMax thresClusterMin offsetCluster facThres precisionDTL useTime
Ts = 0.01; Time = false; %general paras
sigma = 0.01;  winlen=5; thresClusterMax = 1; thresClusterMin = 0.005; offsetCluster = 0; facThres = 2.5; %used in clustering
windowSize = 10; %used in changepoint detection
fixedIntervalLength = 2; precisionDTL = 0.001; useTime = false; %used for DTL
num_var = 1; num_ud = 0; %general paras
max_deriv = 3; %used in changepoint detection
num = 1; x = []; ud = [];

allData = 1:5;
evalData = [4];
%% Changepoint determination and trace setup
tic
for i = allData
    load(['training', int2str(i),'.mat']);
    for j = 1:num_var
        xout(:,j) = 1/max(xout(:,j)) * xout(:,j);
    end
    for deriv = 1:max_deriv
        for curr_var = 1:num_var
            pos_last_deriv = (deriv-1)*num_var + curr_var;
            xout = [xout, [zeros(deriv,1) ; 1/Ts*diff(xout((deriv):end,pos_last_deriv))]];
        end
    end
    %strip info from front bc derivs are not available there
    xout = xout((max_deriv+1):end,:);
    trace_temp = FnDetectChangePoints(xout, num_var);
    trace(num) = trace_temp;
    %all traces are appended (needed for clustering in that form)
    x = [x; trace(num).x];
    ud = [ud; trace(num).ud];
    num = num+1; 
end

%% Determine clustered trace segments

%change num_var so that derivs are also used
num_var = num_var * (1 + 0);
trace = FnClusterSegs(trace, x, ud);
%change back, just in case
num_var = num_var / (1 + 0);

t1 = toc;

%% Create Training Data and run DTL
N = 1;
count = 1;
X = [];
Y = [];
for i = setdiff(allData,evalData)
    [Xnew,Ynew,states] = FnTraceToTrainingData(trace(count),N);
    X = [X; Xnew];
    Y = [Y; Ynew];
    count= count + 1;
end
[Mdl,impure_leaves,num_nodes,learn_time] = FnBuildDT(X,Y);

%% Create Eval Data und evaluate

correct = 0;
false = 0;
Xe = [];
Ye = [];

for i = evalData
    [Xnew,Ynew,states] = FnTraceToTrainingData(trace(i),N);
    Xe = [Xe; Xnew];
    Ye = [Ye; Ynew];
end

for i = 1:size(Xe,1)
    if(predict(Mdl,Xe(i,:)) == Ye(i,1))
        correct = correct + 1;
    else
        false = false + 1;
        disp([mat2str((i+N)*fixedIntervalLength),' in: ',mat2str(Xe(i,:)),' real: ',mat2str(Ye(i,:)),' predict: ',mat2str(predict(Mdl,Xe(i,:)))]);
    end
end

%% ODE Estimation from clustred trace segments
%Limitation: all vars on one segment have the same ODE degree
%            ODE est. has to fail so that lower degree is used

for n =1:length(trace)
    trace(n).labels_trace = [trace(n).labels_trace;0];
end

%change num_var so that derivs are also used
num_var = num_var * (1 + 0);
len_labels = length(trace(1).labels_num);
for label = 1:len_labels
    %todo not only one try but more
    try
        ode(label) = FnEstODE(trace,label);
    catch
        num_var = num_var / (1 + 0);
        num_var = num_var * (1 + 0);
        ode(label) = FnEstODE(trace,label);
        num_var = num_var / (1 + 0);
        num_var = num_var * (1 + 0);
    end
end
%change back, just in case
num_var = num_var / (1 + 0);
t2 = toc;

%% Prediction

offsetPred = 0+1;
sim_x = trace(1).x(offsetPred,:);
lastSwitch = 1;
sim_state = 1;
curr_state = 1;
for j = 1:size(trace(1).chpoints,1)
    if trace(1).chpoints(j,1) > offsetPred
        break;
    end
    lastSwitch = trace(1).chpoints(j,1);
    sim_state = trace(1).labels_trace(j,1);
    curr_state = trace(1).labels_trace(j,1);
end
for i = ((offsetPred-1)/fixedIntervalLength+1):floor(size(trace(1).x)/fixedIntervalLength)
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

hold on;
subplot(1,1,1);
ylim([0 1]);
plot(offsetPred+fixedIntervalLength*(0:(size(sim_x,1)-1)),sim_x(:,1),1:size(trace(1).x,1),trace(1).x(:,1));
hold off;