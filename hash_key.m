function keys = hash_key(feat, tol)
% feat in degrees, tol in degrees

base = round(feat / tol);

keys = {};
idx = 1;

for dx = -1:1
    for dy = -1:1
        for dz = -1:1
            neighbor = base + [dx dy dz];
            keys{idx} = mat2str(neighbor);
            idx = idx + 1;
        end
    end
end
end