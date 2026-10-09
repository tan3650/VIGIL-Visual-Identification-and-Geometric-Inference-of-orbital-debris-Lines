function [R_best, bestMatch, v_rotated] = attitude_determination(v_detected, v_catalog, matches)

if isempty(matches)
    error('No valid matches found.');
end

bestScore = -1; bestErr = inf;
R_best = eye(3);
bestMatch = struct('det_indices',[], 'cat_indices',[]);
v_rotated = [];

for k = 1:length(matches)

    m = matches{k};
    v1 = v_detected(m.det_indices,:);
    v2 = v_catalog(m.cat_indices,:);

    if size(v1,1) < 3, continue; end

    % Kabsch
    H = v1' * v2;
    [U,~,V] = svd(H);
    R = V * U';
    if det(R) < 0, V(:,end) = -V(:,end); R = V * U'; end

    v_all_rot = (R * v_detected')';
    dist = pdist2(v_all_rot, v_catalog);

    score = sum(min(dist,[],2) < 0.01);
    mean_err = mean(min(dist,[],2));

    if score > bestScore || (score == bestScore && mean_err < bestErr)
        bestScore = score;
        bestErr = mean_err;
        R_best = R;
        bestMatch = m;
        v_rotated = v_all_rot;
    end
end

if bestScore < 1
    error('No consistent solution found.');
end

end