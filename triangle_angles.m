function feat = triangle_angles(v)
% v = 3x3 unit vectors

d12 = acosd(dot(v(1,:), v(2,:)));
d23 = acosd(dot(v(2,:), v(3,:)));
d31 = acosd(dot(v(3,:), v(1,:)));

feat = sort([d12, d23, d31]); % SORT = rotation invariant
end