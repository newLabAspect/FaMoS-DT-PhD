function [t_seg,t_cluster,t_charac,t_extract, mse, omega_seg, omega_group] = traceMain(allData,evalData,folder)
    %only global vars needed for this top level program are listed here, if 
    %needed in subfunctions they are listed only there for simplicity
    global num_var num_ud useLMIrefine methodCluster methodTraining offsetCluster
    
    num = 1; x = []; ud = []; xs = [];
    
    addpath(folder);
    % Changepoint Detection only done by proposed algorithm
    addpath(['../../../../ProposedAlgorithm', filesep, 'src']);
    addpath(['../../../../HAutLearn', filesep, 'src']);
    
    %% Changepoint determination and trace setup
    normalization = zeros(num_var+num_ud,1);
    for i = allData
        load(['training', int2str(i),'.mat']);
        for j = 1:num_var
            normalization(j,1) = max(normalization(j,1),max(xout(:,j)));
        end
        for j = 1:num_ud
            normalization(num_var+j,1) = max(normalization(num_var+j,1),max(abs(udout(:,j))));
        end
    end
    num_states = 0;
    tic
    for i = allData
        udout = []; % Needed if system with no input variables present
        load(['training', int2str(i),'.mat']);
        [xout, udout, xout_shifts] = FnShiftAndDiff(xout, udout, normalization);
        trace_temp = FnDetectChangePoints(xout, udout, xout_shifts);
        trace_temp.true_states = states;
        num_states = max([num_states; states]);
        trace_temp.true_chps = chpoints;
        trace(num) = trace_temp;
        %all traces are appended (needed for clustering in that form)
        x = [x; trace(num).x];
        xs = [xs; trace(num).xs];
        ud = [ud; trace(num).ud];
        num = num+1; 
    end
    t_seg = toc;

    trainData = setdiff(allData,evalData);

    % Eval chpoints
    omega_seg = 0.0;
    for i = trainData
        omega_seg = omega_seg + FnEvalChangePoints(trace(i).chpoints, trace(i).true_chps);
    end
    
    %% Determine clustered trace segments
    
    tic;
    
    % Choose which algorithm to use for clustering
    if(methodCluster == 0) % Use DTW for clustering
        useLMIrefine = 0;
        trace = FnClusterSegsFast(trace, x, ud, xs);
    elseif(methodCluster == 1) % Use DTW refined by LMIs for clustering
        useLMIrefine = 1;
        trace = FnClusterSegsFast(trace, x, ud, xs);
    else % Use LMIs for clustering
        trace = FnClusterSegs(trace, xs, ud);
    end
    
    t_cluster = toc;

    % Remove changepoints that are not associated with a system mode switch
    trace = FnCleanChangePoints(trace);
    
    save("trace_before_ode", "trace");
    tic
    ode = FnEstODE(trace(trainData));
    t_charac = toc;

    %Cluster and Characterization Eval
    mse_vec = zeros(length(trainData),num_var*(1+offsetCluster));
    mse_per_group = [];
    for i = 1:length(trainData)
        [mse_vec(i,:), ~, mpg, gs] = FnEvalFlowAccuracy(trace(trainData(i)), ode);
        for g = 1:length(mpg)
            if length(mse_per_group) < g
                mse_per_group(g) = mpg(g);
                group_sizes(g) =  gs(g);
            else
                mse_per_group(g) = mse_per_group(g) + mpg(g);
                group_sizes(g) = group_sizes(g) + gs(g);
            end
        end
    end
    omega_group = sum(mse_per_group .* group_sizes) / sum(group_sizes) * (1 + abs(num_states - length(mse_per_group)));
    
    %% Training and Eval
    
    % Choose which algorithm to use for training
    if(methodTraining == 0) % Use DTL for training
        
        %Extraction
        X = [];
        Y = [];
        for i = trainData
            [Xnew,Ynew,~] = FnTraceToTrainingData(trace(i));
            X = [X; Xnew];
            Y = [Y; Ynew];
        end

        tic
        [Mdl,~,~,~] = FnBuildDT(X,Y);
        t_extract = toc;
    
        % Actual Eval
        mse_vec = zeros(length(evalData),num_var*(1+offsetCluster));
        for i = 1:length(evalData)
            gt_trace = trace(evalData(i));
            %Prediction
            [sim_trace] = FnPredictTraceDTL(gt_trace,Mdl,ode);
            for j = 1:size(mse_vec,2)
                mse_vec(i,j) = mean((sim_trace.x(:,j)-gt_trace.xs(:,j)).^2);
            end
            traces(i) = sim_trace;
        end
        gt_traces = trace(evalData);
        mse = mean(mse_vec,1);
        save(strcat("model"), "Mdl", "ode", "traces", "gt_traces", "mse_vec");

    else % Use PTA for training
        global eta lambda gamma
        
        if size(ode, 2) > 1
            for n =1:length(trace)
                trace(n).labels_trace = [trace(n).labels_trace;0];
            end
            
            tic;
            % LI Estimation given parameters
            [trace_train,label_guard] = FnLI(trace(trainData), eta, lambda, gamma);
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
            t_extract = toc;
    
            xmlstruct = readstruct([folder, filesep, 'automata_learning.xml']);
            
            conditions = FnParseHA(xmlstruct,ode,label_guard);
        else
            model = [];
            conditions = [];
            t_extract = NaN;
        end

        mse_vec = zeros(length(evalData),num_var*(1+offsetCluster));
        for i = 1:length(evalData)
            gt_trace = trace(evalData(i));
            %Prediction
            [sim_trace] = FnPredictTraceHA(gt_trace,conditions,ode);
            for j = 1:size(mse_vec,2)
                mse_vec(i,j) = mean((sim_trace.x(:,j)-gt_trace.xs(:,j)).^2);
            end
            traces(i) = sim_trace;
        end
        gt_traces = trace(evalData);
        mse = mean(mse_vec,1);
        save("model", "conditions", "ode", "traces", "gt_traces", "mse_vec");
    end
end
