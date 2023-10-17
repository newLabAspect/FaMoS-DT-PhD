function [xout, udout, xout_shifts] = FnShiftAndDiff(xout, udout, norm_coeff)
% FnNormAndDiff duplicates and shifts output variables for easier discrete
% state space calculations and calculates deriatives up to max_deriv
    global num_ud num_var max_deriv Ts

    % Norm xout
    for j = 1:num_var
        xout(:,j) = 1/norm_coeff(j,1) * xout(:,j);
    end

    % Calculate shifted duplicates
    xout_shifts = zeros(size(xout,1),num_var*(1+max_deriv));
    for shift = 0:max_deriv
        xout_shifts(:,(1+shift*num_var):((shift+1)*num_var)) = [zeros(shift,num_var); xout(1:(end-shift),1:num_var)];
    end

    % Calculate derivatives up to selected degree, structure of xout is
    % [x_1,...,x_n,x_1',...,x_n',...] thus curr_var in inner loop
    for deriv = 1:max_deriv
        for curr_var = 1:num_var
            pos_last_deriv = (deriv-1)*num_var + curr_var;
            xout = [xout, [zeros(deriv,1) ; 1/Ts*diff(xout((deriv):end,pos_last_deriv))]];
        end
    end

    % Strip entries from front of outputs bc derivs/shifts are not available there
    xout = xout((max_deriv+1):end,:);
    xout_shifts = xout_shifts((max_deriv+1):end,:);

    % Normalize inputs using normalization factors (Derivs not needed)
    for j = 1:num_ud
        udout(:,j) = 1/norm_coeff(num_var+j,1) * udout(:,j);
    end

    % Strip entries from front of inputs to match with output vars
    if num_ud ~= 0
        udout = udout((max_deriv+1):end,:);
    end
end