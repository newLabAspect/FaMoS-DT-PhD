%% Intital Setup and parameter settings
clc
clear

addpath(['InterOperability']);

%% General paras
global num_var num_ud Ts Time methodCluster methodTraining windowSize max_deriv offsetCluster
num_var = 1; num_ud = 0;
methodCluster = 0; % 0: DTW, 1: DTW & LMI, 2: LMI
methodTraining = 0; % 0: DTL, 1: PTA
Ts = 0.01; Time = false;
% Changepoint detection paras
windowSize = 10; max_deriv = 3;
% Up to which shift to consider in LMI, Automata Learning
% On which derivative to cluster using DTW comparisons
offsetCluster = 1;

%% Clustering Paras
% LMI paras (can comment out if only DTW used)
global sigma winlen
sigma = 0.0000001;  winlen=1; % LMI paras
% DTW paras (can comment out if only LMIs used)
global thresClusterMax thresClusterMin facThres
thresClusterMax = 0.1; thresClusterMin = 0.05; facThres = 7.5;
% Slightly Off Paras:
%thresClusterMax = 0.04; thresClusterMin = 0; facThres = 2.5;

%% Training Paras
% PTA paras (can comment out if DTL is used)
global eta lambda gamma tolLI
eta = 100000; % number of iterations 
lambda = 0.1; % tolerance 
gamma = 5; %the least number of inlayers
tolLI = 0.0037; %tolerance in evaluation of LIs
% DTL paras (can comment out if PTA is used)
global fixedIntervalLength precisionDTL useTime
fixedIntervalLength = 1; precisionDTL = 0.0001; useTime = false;

%% Vary Paras over time
global variedMetric variedMetricSteps
variedMetric = 4; % -1: No parameter is varied
%DTL: 4: trainingSetSize
%PTA: 0: eta 1: lambda 2: gamma 3: toLi 4: trainingSetSize
variedMetricSteps = linspace(0.25,0.25,1);

%% Actual execution
allData = 1:10;
evalData = [2,8];

[correct,false,t_cluster,t_train,trace,ClusterCorrect,ClusterFalse,pred_trace,confDeg] = traceMain(allData,evalData,['ExampleSystems', filesep, 'TwoStateHybridAutomaton']);