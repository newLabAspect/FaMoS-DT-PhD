%% Intital Setup and parameter settings
clc
clear

addpath(['ExampleSystems', filesep, 'MultiVariablesHeatingSystem']);
addpath(['HAutLearn', filesep, 'src']);
global sigma num_var num_ud winlen Ts Time
Ts = 0.01; Time = false;
sigma = 0.00005;  winlen=5;
num_var = 2; num_ud = 0;
num = 1; x = []; ud = [];

allData = 1:5;
evalData = [4,5];
%% Changepoint determination and trace setup
tic
for i = allData
    load(['training', int2str(i),'.mat']);
    trace_temp = FnProcessNoiseData(xout, num_var);
    trace(num) = trace_temp;
    x = [x; trace(num).x];
    ud = [ud; trace(num).ud];
    num = num+1; 
end

%% Parameter settings for edge conditioning
t0 = toc;
eta = 100000; % number of iterations 
lambda = 0.1; % tolerance 
gamma = 10; %the least number of inlayers

%% Determine clustered trace segments

%trace = FnClusterSegs(trace, x, ud);

t1 = toc;

%% ODE Estimation from clustred trace segments

for n =1:length(trace)
    trace(n).labels_trace = [trace(n).labels_trace;0];
end

ode = FnEstODE(trace(setdiff(allData,evalData)));
t2 = toc;

%% LI Estimation given parameters
[trace_train,label_guard] = FnLI(trace(setdiff(allData,evalData)), eta, lambda, gamma);
t3 = toc;

%% Setup PTA given LIs and ODEs
pta_trace = FnPTA(trace_train);
pta_trace = pta_filter(pta_trace);
t4 = toc;

%% Generate Final Automaton model

FnGenerateHyst(['ExampleSystems', filesep, 'MultiVariablesHeatingSystem', filesep, 'automata_learning'],label_guard, num_var, ode, pta_trace);

xmlstruct = readstruct(['ExampleSystems', filesep, 'MultiVariablesHeatingSystem', filesep, 'automata_learning.xml']);

[correct,false] = FnEvaluate(trace(evalData),xmlstruct,ode,0.001);

t5 = toc;


%% Additional functions (case dependent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pta_trace_new = pta_filter(pta_trace)
    % remove false pta
    label1s = extractfield(pta_trace,'label1');
    label2s = extractfield(pta_trace,'label2');
    id1s = extractfield(pta_trace,'id1');
    id2s = extractfield(pta_trace,'id2');
    nn = 1;
    while true
        if nn>length(pta_trace)
            break;
        end
        %flag1 = pta_trace(nn).times<=2;
        flag2 = ~ismember(pta_trace(nn).label1, label2s);
        flag3 = ~ismember(pta_trace(nn).label2, label1s);
        flag4 = ~ismember(pta_trace(nn).id1, id1s);
        flag5 = ~ismember(pta_trace(nn).id2, id2s);
        
        if (~ismember(pta_trace(nn).id1, id2s)) && (pta_trace(nn).id1~=1) || ~ismember(pta_trace(nn).id2, id1s)
            pta_trace(nn) = [];
        else
            nn = nn+1;
        end
    end
    pta_trace_new = pta_trace;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trace = FnProcessNoiseData(xout, num_var)

    chpoints = [];
    for i = 1:num_var
        chpoints = union(chpoints, changepoint(xout(:,i)));
    end
    
    % remove redundant chpoints
    chpoints = filterindx(chpoints);
    xout_reduced= xout(:, 1:num_var);

    trace.x = xout_reduced;
    trace.chpoints = chpoints;
    trace.ud = [];
    trace.labels_num = []; 
    trace.labels_trace = [];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indx = changepoint(values)
    len = 2;
    diffs = diff(values(1:end-len,:)-values((len+1):end,:));
    indx = find(diffs<=-0.18|diffs>=0.25)+len;
    indx = union(1,[indx; length(values)]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indx = filterindx(indx)
    n = 1;
    windw = 10; 
    while true
        if n >= length(indx)
            break;
        end
        id1 = indx(n);
        while true
            if n+1 >= length(indx)
                break;
            end
            id2 = indx(n+1);
            if id2-id1<=windw
                indx(n+1) = [];
            else
                break;
            end
        end
        n = n+1;
    end
end