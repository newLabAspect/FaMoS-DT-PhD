function rounded = FnRoundToInterval(vals,step)
% FnRoundToInterval rounds all entries in the input to the nearest
% quantization levels. All levels are placed equidistantly around zero.
    
    rounded = zeros(size(vals));
    for i = 1:size(vals,1) % Row
        for j = 1:size(vals,2) % Column
            val = vals(i,j);
            % Transform value onto integer grid
            val = val * 1/step;
            % Round on integer grod
            val = round(val,0);
            % Transform back
            rounded(i,j) = val * step;
        end
    end
end