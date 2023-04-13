%% Intital Setup and parameter settings
clc
clear

addpath(['ExampleSystems', filesep, 'ComplexTankSystem']);
addpath(['ProposedAlgorithm', filesep, 'src']);
global sigma num_var num_ud winlen Ts Time windowSize fixedIntervalLength max_deriv thresClusterMax thresClusterMin offsetCluster facThres precisionDTL useTime useLMIrefine
Ts = 0.01; Time = false; %general paras
sigma = 0.000001;  winlen=1; thresClusterMax = 1; thresClusterMin = 0.01; offsetCluster = 0; facThres = 2.5; useLMIrefine = 0;%used in clustering
Nsteps = 1; rangePrecision = 0.0001; %linspace(0.01,0.0001,500); %used in eval
windowSize = 10; %used in changepoint detection
fixedIntervalLength = 1; precisionDTL = 0.001; useTime = true; %used for DTL

num_var = 3; num_ud = 0; %general paras
max_deriv = 3; %used in changepoint detection
num = 1; x = []; ud = [];

allData = 1:10;
evalData = [2,6];
%% Changepoint determination and trace setup
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

tic;
%change num_var so that derivs are also used
num_var = num_var * (1 + 0);
trace = FnClusterSegs(trace, x, ud);
%change back, just in case
num_var = num_var / (1 + 0);

t_cluster = toc;

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
        
        tic;
        [Mdl,impure_leaves,num_nodes,learn_time] = FnBuildDT(X(1:k,:),Y(1:k,:));
        t_train = toc;

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
    writematrix([rangePrecision',falseAll],['ExampleSystems', filesep, 'ComplexTankSystem', filesep, 'COTS_quant_time.csv']);
elseif (Nsteps ~= 1)
    writematrix([round(linspace(1,size(X,1),Nsteps),0)',falseAll],['ExampleSystems', filesep, 'ComplexTankSystem', filesep, 'COTS_size_time_smallquant.csv']);
end