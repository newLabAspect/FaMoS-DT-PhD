function rounded = FnRoundToInterval(vals,step)
    rounded = zeros(size(vals));
    for i = 1:size(vals,1) %row
        for j = 1:size(vals,2) %column
            val = vals(i,j);
            val = val * 1/step;
            val = round(val,0);
            rounded(i,j) = val * step;
        end
    end
end