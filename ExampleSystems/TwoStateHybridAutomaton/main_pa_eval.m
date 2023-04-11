%% Intital Setup and parameter settings
clc
clear

addpath(['ExampleSystems', filesep, 'TwoStateHybridAutomaton']);
addpath(['ProposedAlgorithm', filesep, 'src']);
global sigma num_var num_ud winlen Ts Time windowSize fixedIntervalLength max_deriv thresClusterMin thresClusterMax offsetCluster facThres precisionDTL useTime
Ts = 0.01; Time = false; %general paras
sigma = 0.000005;  winlen=5; thresClusterMin = 0; thresClusterMax = 0.04; offsetCluster = 1; facThres = 2.5; %used in clustering
Nsteps = 1; rangePrecision = 0.0001; %linspace(0.01,0.0001,500); %used in eval
windowSize = 10; %used in changepoint detection
fixedIntervalLength = 1; precisionDTL = 0.001; useTime = true; %used for DTL

num_var = 1; num_ud = 0; %general paras
max_deriv = 3; %used in changepoint detection
num = 1; x = []; ud = [];

allData = 1:10;
evalData = [2,8];
%% Changepoint determination and trace setup
tic
normalization = zeros(num_var,1); % obviously not ideal but works
for i = allData
    load(['training', int2str(i),'.mat']);
    for j = 1:num_var
        normalization(j,1) = max(normalization(j,1),max(xout(:,j)));
    end
end
for i = allData
    load(['training', int2str(i),'.mat']);
    for j = 1:num_var
        xout(:,j) = 1/normalization(j,1) * xout(:,j);
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

%% Manual fixup (in future done by LMI)

for i = 1:10
    for j = 1:length(trace(i).labels_trace)
        if(trace(i).labels_trace(j) == 7)
            trace(i).labels_trace(j) = 1;
        elseif(trace(i).labels_trace(j) == 13)
            trace(i).labels_trace(j) = 1;
        elseif(trace(i).labels_trace(j) == 26)
            trace(i).labels_trace(j) = 1;
        end
    end
end

%% Eval

%legacy eval paras
N = 1; % how many past states should be used to predict

%To save results
correctAll = [];
falseAll = [];
for h = 1:length(rangePrecision)
    precisionDTL = rangePrecision(1,h);

    %compute eval data
    Xe = [];
    Ye = [];
    
    for i = evalData
        [Xnew,Ynew,states] = FnTraceToTrainingData(trace(i),N);
        Xe = [Xe; Xnew];
        Ye = [Ye; Ynew];
    end

    %compute maximal training data
    X = [];
    Y = [];
    count = 1;
    for i = setdiff(allData,evalData)
        [Xnew,Ynew,states] = FnTraceToTrainingData(trace(count),N);
        X = [X; Xnew];
        Y = [Y; Ynew];
        count= count + 1;
    end

    for k = round(linspace(1,size(X,1),Nsteps),0)
        disp([num2str(round(k/size(X,1)*h/length(rangePrecision)*100,2)),' % done']);
        %DTL
        
        [Mdl,impure_leaves,num_nodes,learn_time] = FnBuildDT(X(1:k,:),Y(1:k,:));
        
        % Actual Eval
        correct = 0;
        false = 0;
        
        for i = 1:size(Xe,1)
            if(predict(Mdl,Xe(i,:)) == Ye(i,1))
                correct = correct + 1;
            else
                false = false + 1;
                %disp([mat2str((i+N)*fixedIntervalLength),' in: ',mat2str(Xe(i,:)),' real: ',mat2str(Ye(i,:)),' predict: ',mat2str(predict(Mdl,Xe(i,:)))]);
            end
        end
        correctAll = [correctAll; correct];
        falseAll = [falseAll; false];
    end
end
if(length(rangePrecision) ~= 1)
    writematrix([rangePrecision',falseAll],['ExampleSystems', filesep, 'TwoStateHybridAutomaton', filesep, '2SHA_quant_time.csv']);
elseif (Nsteps ~= 1)
    writematrix([round(linspace(1,size(X,1),Nsteps),0)',falseAll],['ExampleSystems', filesep, 'TwoStateHybridAutomaton', filesep, '2SHA_size_time_smallquant.csv']);
end