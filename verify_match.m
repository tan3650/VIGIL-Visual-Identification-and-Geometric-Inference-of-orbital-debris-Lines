function [R, err] = verify_match(v_detected, v_catalog, det_idx, cat_idx)
v_det = v_detected(det_idx, :)';
v_cat = v_catalog(cat_idx, :)';

% Wahba using SVD
B = v_cat * v_det';
[U,~,V] = svd(B);
R = U * diag([1 1 det(U*V')]) * V';

% Apply rotation to all detected vectors
v_cat_est = (R * v_detected')';
dist = pdist2(v_cat_est, v_catalog);
err = mean(min(dist, [], 2));  % Mean reprojection error
end