function [correctAll,falseAll,t_cluster,t_train,trace,ClusterCorrect,ClusterFalse] = evalMain(allData,evalData,folder)
    %only global vars needed for this top level program are listed here, if 
    %needed in subfunctions they are listed only there for simplicity
    global num_var num_ud Ts max_deriv useLMIrefine methodCluster methodTraining variedMetricSteps variedMetric offsetCluster
    
    num = 1; x = []; ud = [];
    
    addpath(folder);
    % Changepoint Detection only done by proposed algorithm
    addpath(['ProposedAlgorithm', filesep, 'src']);
    addpath(['HAutLearn', filesep, 'src']);
    
    %% Changepoint determination and trace setup
    
    % Normalize trace to [-1,+1] by computing maximum absolute trace value
    normalization = zeros(num_var+num_ud,1);
    for i = allData
        load(['training', int2str(i),'.mat']);
        for j = 1:num_var
            normalization(j,1) = max(normalization(j,1),max(abs(xout(:,j))));
        end
        for j = 1:num_ud
            normalization(num_var+j,1) = max(normalization(num_var+j,1),max(abs(udout(:,j))));
        end
    end

    % Carry out normalization, compute derivatives and detect changepoints
    % Results are written into trace datastructure (incl. ground truths)
    for i = allData
        udout = []; % Needed if system with no input variables present
        load(['training', int2str(i),'.mat']);
        [xout, udout] = FnNormAndDiff(xout, udout, normalization);
        trace_temp = FnDetectChangePoints(xout, udout);
        % Ground truths for states and changepoints
        trace_temp.true_states = states;
        trace_temp.true_chps = chpoints;
        trace(num) = trace_temp;
        % All trace values are combined into one array as this form is
        % needed in the clustering stage
        x = [x; trace(num).x];
        ud = [ud; trace(num).ud];
        num = num+1; 
    end
    
    %% Determine clustered trace segments
    
    tic;

    % Choose which algorithm to use for clustering
    if(methodCluster == 0) % Use DTW for clustering
        useLMIrefine = 0;
        trace = FnClusterSegsFast(trace, x, ud);
    elseif(methodCluster == 1) % Use DTW refined by LMIs for clustering
        useLMIrefine = 1;
        trace = FnClusterSegsFast(trace, x, ud);
    else % Use LMIs for clustering
        trace = FnClusterSegs(trace, x, ud);
    end
    
    t_cluster = toc;
    
    % Remove changepoints that are not associated with an system mode switch
    trace = FnCleanChangePoints(trace);

    %Eval clusters
    ClusterCorrect = 0;
    ClusterFalse = 0;
    for i = 1:length(trace)
        [cTemp, fTemp] = FnEvalCluster(trace(i).labels_trace,trace(i).true_states,trace(i).true_chps);
        ClusterCorrect = ClusterCorrect + cTemp;
        ClusterFalse = ClusterFalse + fTemp;
    end
    
    %% Training and short Eval
    
    % Choose which algorithm to use for training
    if(methodTraining == 0) % Use DTL for training
        global precisionDTL
        %legacy eval paras
        N = 1; % how many past states should be used to predict
        
        correctAll = [];
        falseAll = [];
        t_train = [];

        for j = 1:length(variedMetricSteps)
            %vary selected metric
            if(variedMetric == 0)
                precisionDTL = variedMetricSteps(1,j);
            end

            %compute eval data
            Xe = [];
            Ye = [];
            
            for i = evalData
                [Xnew,Ynew,states] = FnTraceToTrainingData(trace(i));
                Xe = [Xe; Xnew];
                Ye = [Ye; Ynew];
            end
            
            tic;
            %compute maximal training data
            X = [];
            Y = [];
            count = 1;
            for i = setdiff(allData,evalData)
                [Xnew,Ynew,states] = FnTraceToTrainingData(trace(count));
                X = [X; Xnew];
                Y = [Y; Ynew];
                count= count + 1;
            end
            
            %vary selected metric
            if(variedMetric == 1)
                X = X(1:round(variedMetricSteps(1,j)*size(X,1),0),:);
                Y = Y(1:round(variedMetricSteps(1,j)*size(Y,1),0),:);
            end

            [Mdl,impure_leaves,num_nodes,learn_time] = FnBuildDT(X,Y);
            t_train = [t_train; toc];
        
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
            if(variedMetric == -1)
                break;
            end
            disp([num2str(round(j/length(variedMetricSteps)*100,2)),' % done']);
        end
    else % Use PTA for training
        global eta lambda gamma tolLI

        for n =1:length(trace)
            trace(n).labels_trace = [trace(n).labels_trace;0];
        end

        correctAll = [];
        falseAll = [];
        t_train = [];

        trace_train_short = trace(setdiff(allData,evalData));
        max_len_trace = 0;
        for p = 1:length(trace_train_short)
            max_len_trace = max_len_trace + size(trace_train_short(p).x,1);
        end
        
        %maybe move for higher throughput
        for j = 1:length(variedMetricSteps)
            %vary selected metric
            if(variedMetric == 0)
                eta = variedMetricSteps(1,j);
            elseif(variedMetric == 1)
                lambda = variedMetricSteps(1,j);
            elseif(variedMetric == 2)
                gamma = variedMetricSteps(1,j);
            elseif(variedMetric == 3)
                tolLI = variedMetricSteps(1,j);
            elseif(variedMetric == 4)
                used_len = round(variedMetricSteps(1,j) * max_len_trace,0);
                trace_train_short = trace(setdiff(allData,evalData));
                labels_num = [];
                for p = 1:length(trace_train_short)
                    if(used_len == 0)
                        % all remaining trace entries need to be deleted
                        trace_train_short(p:length(trace_train_short)) = [];
                        break;
                    elseif(used_len < length(trace_train_short(p).x))
                        % trim x, ud entries out of trace
                        trace_train_short(p).x = trace_train_short(p).x(1:used_len,:);
                        if num_ud ~= 0
                            trace_train_short(p).ud = trace_train_short(p).ud(1:used_len,:);
                        end
                        % trim changepoints out of trace
                        toDelete = find(trace_train_short(p).chpoints >= used_len);
                        trace_train_short(p).chpoints(toDelete) = [];
                        % trim cluster ids out of trace
                        trace_train_short(p).labels_trace = trace_train_short(p).labels_trace(1:length(trace_train_short(p).chpoints));
                        % add ending chpoint to trace
                        trace_train_short(p).chpoints = [trace_train_short(p).chpoints; used_len];
                        % add terminating 0 to cluster ids
                        trace_train_short(p).labels_trace = [trace_train_short(p).labels_trace;0];
                        % trace data structure is completed
                        used_len = 0;
                    else
                        % trace data structure still needs entries to achieve target size
                        used_len = used_len - length(trace_train_short(p).x);
                    end
                    %to exclude to trailing 0
                    labels_num = unique([labels_num; trace_train_short(p).labels_trace(1:(end-1))]);
                end
                for p = 1:length(trace_train_short)
                    trace_train_short(p).labels_num = labels_num;
                end
            end

            if(j == 1 || variedMetric ~= 3)
                %ODE estimation excluded from time, because both do them
                ode = FnEstODE(trace_train_short);
                
                tic;
                % LI Estimation given parameters
                [trace_train,label_guard] = FnLI(trace_train_short, eta, lambda, gamma);
                % Extend label_guard by zeros for considered derivatives so
                % that generated automaton xml file is clean (currently no
                % transitions based on derivatives are allowed)
                for k = 1:length(label_guard)
                    curr_label_guard = cell2mat(label_guard(k));
                    label_guard(k) = {[curr_label_guard(1:num_var); zeros(offsetCluster*num_var,1); curr_label_guard((num_var+1):end)]};
                end
    
                % Setup PTA given LIs and ODEs
                pta_trace = FnPTA(trace_train);
                % Can solve bugs if you comment out next line. But why ?
                %pta_trace = pta_filter(pta_trace);
                
                % Generate Final Automaton model
                
                FnGenerateHyst([folder, filesep, 'automata_learning'],label_guard, num_var*(1+offsetCluster), num_ud, ode, pta_trace);
                t_train = [t_train; toc]; % How to measure this time should be discussed
            end % What happens with variables in scope? Seems to work...

            xmlstruct = readstruct([folder, filesep, 'automata_learning.xml']);
            
            conditions = FnParseHA(xmlstruct,ode,label_guard);
            [correct,false] = FnEvaluateHA(trace(evalData),conditions,tolLI);
            correctAll = [correctAll; correct];
            falseAll = [falseAll; false];
            if(variedMetric == -1)
                break;
            end
            disp([num2str(round(j/length(variedMetricSteps)*100,2)),' % done']);
        end
    end
end

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
