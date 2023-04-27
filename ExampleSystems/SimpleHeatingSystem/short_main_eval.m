%% Intital Setup and parameter settings
clc
clear

addpath(['InterOperability']);

%% General paras
global num_var num_ud Ts Time methodCluster methodTraining windowSize max_deriv
num_var = 1; num_ud = 0;
methodCluster = 0; % 0: DTW, 1: DTW & LMI, 2: LMI
methodTraining = 1; % 0: DTL, 1: PTA
Ts = 0.01; Time = false;
% Changepoint detection paras
windowSize = 10; max_deriv = 3;

%% Clustering Paras
% LMI paras (can comment out if only DTW used)
global sigma winlen
sigma = 0.00005;  winlen=5; % LMI paras
% DTW paras (can comment out if only LMIs used)
global thresClusterMax thresClusterMin offsetCluster facThres
thresClusterMax = 0.1; thresClusterMin = 0.01; offsetCluster = 0; facThres = 2.5;

%% Training Paras
% PTA paras (can comment out if DTL is used)
global eta lambda gamma tolLI
eta = 100000; % number of iterations 
lambda = 0.1; % tolerance 
gamma = 10; %the least number of inlayers
tolLI = 0.0022; %tolerance in evaluation of LIs
% DTL paras (can comment out if PTA is used)
global fixedIntervalLength precisionDTL useTime
fixedIntervalLength = 1; precisionDTL = 0.0001; useTime = true;

%% Vary Paras over time
global variedMetric variedMetricSteps
variedMetric = -1; % -1: No parameter is varied
%DTL: 0: precisionDTL 1: trainingSetSize
%PTA: 0: eta 1: lambda 2: gamma 3: toLi
variedMetricSteps = linspace(0.01,0.001,100);

%% Actual execution
allData = 1:10;
evalData = [2,5];

[correct,false,t_cluster,t_train,trace] = evalMain(allData,evalData,['ExampleSystems', filesep, 'SimpleHeatingSystem']);