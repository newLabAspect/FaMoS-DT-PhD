function [confDeg] = FnConformanceTrace(trace, est_trace)
% FnConformanceTrace computes the (T,J,tau_max,eps) closeness between two
% traces (usually preexisting and simulated)
%   This function is derived from accuracy_evaluation.m of HAutLearn 
    global num_var Ts

    % Compute tau_max and bring changepoint set to same length
    est_chp = est_trace.chpoints;
    chp = trace.chpoints;
    len_chp = min(length(est_chp), length(chp));
    
    est_chp2 = est_chp(1:len_chp);
    chp2 = chp(1:len_chp);
    
    tau_max = max(abs(est_chp2-chp2));

    % Compute eps
    for j = 2:len_chp
        startj = chp(j-1);
        endj = chp(j);
        startj_est = est_chp(j-1);
        endj_est = est_chp(j);
        
        % Extract currently considered segment
        xs = trace.x(startj: endj,1:num_var);
        xs_est = est_trace.x(startj_est: endj_est,1:num_var);
        
        % Mitigate ordering effects by executing in both orders
        eps = errorAlignedPoints(xs,xs_est,tau_max);
        eps = max(eps, errorAlignedPoints(xs_est,xs,tau_max));
    end

    % Return Conformance Degree 
    confDeg = [length(trace.x)*Ts, len_chp, tau_max*Ts, eps];
end

function [eps] = errorAlignedPoints(xs,xs_est,tau_max)
% errorAlignedPoints computes the maximal distance between matching points
% from two traces allowing from up to tau_max as a timeshift
    global num_var

    % Variable to track error
    eps = 0.0;

    % Go over all points from original trace
    for h = 1:size(xs,1)
        % Due to possible timeshift these indices need to be considered for
        % possible matching points from estimated trace
        start_temp = h-tau_max;
        end_temp = h+tau_max;
        if start_temp<1
            start_temp = 1;
        end
        if end_temp>size(xs_est,1)
            end_temp = size(xs_est,1);
        end

        % Compute errors for each matching and choose minimal one
        xs_temp = xs(h,1:num_var)-xs_est(start_temp:end_temp,1:num_var);
        e_temp = min(diag(xs_temp*xs_temp'));
        eps = max(eps, sqrt(e_temp));
    end
end