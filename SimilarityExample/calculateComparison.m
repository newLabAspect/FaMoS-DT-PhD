hold off;
trace = readmatrix("SimilarityExample\seg_data.csv"); % t, y, deviation, y2, y3, dy, deviation_derivative

plot(trace(:,2))
hold on;
plot(trace(:,4))
plot(trace(:,5))

seg_1 = trace(1:100,2);
seg_2 = trace(100:600,2);
seg_3 = trace(600:700,2);
seg_4 = trace(700:800,2);
seg_5 = trace(800:end,2);

computeSimilarity(seg_1,seg_2)
computeSimilarity(seg_1,seg_3)
computeSimilarity(seg_1,seg_4)
computeSimilarity(seg_1,seg_5)

computeSimilarity(seg_2, seg_3)
computeSimilarity(seg_2, seg_4)
computeSimilarity(seg_2, seg_5)

computeSimilarity(seg_3, seg_4)
computeSimilarity(seg_3, seg_5)

computeSimilarity(seg_4, seg_5)

function sim = computeComparison(seg_1,seg_2)
% computeComparison returns the similarity index of two segments
%   The similarity index is computed using a DTW comparison and combining
%   the distance and diagonality metric into a single one
    [dist, i_x, i_y] = dtw(seg_1,seg_2);
    diag = corrcoef(i_x,i_y);
    diag = diag(1,2);
    dist = dist/size(i_x,2);
    % Diag can be NaN if one signal is a constant resulting in
    % zero variance, in this case no similarity is given
    if isnan(diag)
        diag = 0.0;
    end
    % Combine both metrics into a single one
    sim = 0.5*dist+0.5*(1-diag);
end

function combined_metric = computeSimilarity(seg_1,seg_2)
    common_len = min(size(seg_1,1),size(seg_2,1));
    
    % Compute DTW-Comparison for shortened, at END aligned signals
    sim_index_end = computeComparison(seg_1((end-common_len+1):end,1)-seg_1(end,1),seg_2((end-common_len+1):end,1)-seg_2(end,1));
    
    % Compute DTW-Comparison for shortened, at START aligned signals
    sim_index_start = computeComparison(seg_1(1:common_len,1)-seg_1(1,1),seg_2(1:common_len,1)-seg_2(1,1));

    % Select the better similarity index and save it in
    % symmetrical extended matrix
    if (sim_index_start <= sim_index_end)
        combined_metric = sim_index_start;
    else
        combined_metric = sim_index_end;
    end

end