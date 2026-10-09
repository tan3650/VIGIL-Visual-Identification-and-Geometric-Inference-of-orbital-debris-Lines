% --- MODIFIED triangle_features.m with angles ---
function feat = triangle_features(v)
a = norm(v(2,:) - v(1,:));
b = norm(v(3,:) - v(2,:));
c = norm(v(1,:) - v(3,:));
s = (a + b + c)/2;
area = sqrt(max(s*(s-a)*(s-b)*(s-c), 0));

if area == 0
    feat = [0, 0, 0, 0, 0, 0];
else
    % Normalize side lengths by 2*area (original method)
    lengths = sort([a, b, c] / (2 * area));

    % Compute angles using Law of Cosines
    alpha = acosd((b^2 + c^2 - a^2)/(2 * b * c));
    beta  = acosd((a^2 + c^2 - b^2)/(2 * a * c));
    gamma = 180 - alpha - beta;
    angles = sort([alpha, beta, gamma] / 180);  % normalize to [0, 1]

    % Combine length-based and angle-based features
    feat = [lengths, angles];
end
end