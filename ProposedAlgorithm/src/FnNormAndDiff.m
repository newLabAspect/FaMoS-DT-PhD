function xout = FnNormAndDiff(xout, norm_coeff)
% FnNormAndDiff scales the output variables using the normalization
% coefficients and calculates deriatives up to max_deriv
    global num_var max_deriv Ts

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
end