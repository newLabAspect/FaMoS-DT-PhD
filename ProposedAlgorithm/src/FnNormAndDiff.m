function [xout, udout] = FnNormAndDiff(xout, udout, norm_coeff)
% FnNormAndDiff scales the output variables using the normalization
% coefficients and calculates deriatives up to max_deriv
    global num_ud num_var max_deriv Ts

    % Normalize using normalization factor
    for j = 1:num_var
        xout(:,j) = 1/norm_coeff(j,1) * xout(:,j);
    end

    % Calculate derivatives up to selected degree, structure of xout is
    % [x_1,...,x_n,x_1',...,x_n',...] thus curr_var in inner loop
    for deriv = 1:max_deriv
        for curr_var = 1:num_var
            pos_last_deriv = (deriv-1)*num_var + curr_var;
            xout = [xout, [zeros(deriv,1) ; 1/Ts*diff(xout((deriv):end,pos_last_deriv))]];
        end
    end

    % Strip entries from front bc derivs are not available there
    xout = xout((max_deriv+1):end,:);

    % Inputs are currently only normalized
    for j = 1:num_ud
        udout(:,j) = 1/norm_coeff(num_var+j,1) * udout(:,j);
    end
    % And shortened to match with output vars
    if num_ud ~= 0
        udout = udout((max_deriv+1):end,:);
    end
end