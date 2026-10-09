function triangle_db = build_triangle_db(v_catalog, tol)

idx = nchoosek(1:size(v_catalog, 1), 3);
num_tri = size(idx, 1);

% Preallocate storage
keys_all = cell(num_tri, 1);
tris_all = cell(num_tri, 1);

% PARALLEL
parfor i = 1:num_tri
    tri = v_catalog(idx(i,:), :);

    feat = triangle_angles(tri);
    keys = hash_key(feat, tol);

    keys_all{i} = keys;
    tris_all{i} = idx(i,:);
end

% SEQUENTIAL MAP BUILD 
triangle_db = containers.Map;

for i = 1:num_tri
    keys = keys_all{i};
    tri_idx = tris_all{i};

    for k = 1:length(keys)
        key = keys{k};

        if isKey(triangle_db, key)
            temp = triangle_db(key);
            temp{end+1} = tri_idx;
            triangle_db(key) = temp;
        else
            triangle_db(key) = {tri_idx};
        end
    end
end
end