function matches = pattern_matching(v_detected, triangle_db, v_catalog, tol)

det_idx = nchoosek(1:size(v_detected,1),3);
matches = {};
count = 0; max_matches = 100;

for i = 1:size(det_idx,1)

    tri = v_detected(det_idx(i,:),:);
    feat = triangle_angles(tri);
    keys = hash_key(feat, tol);

    for k = 1:length(keys)
        key = keys{k};

        if isKey(triangle_db,key)
            cat_tris = triangle_db(key);

            for j = 1:length(cat_tris)

                match.det_indices = det_idx(i,:);
                match.cat_indices = cat_tris{j};

                [R, err] = verify_match(v_detected, v_catalog, ...
                    match.det_indices, match.cat_indices);

                if err < 0.02

                    det_extra = setdiff(1:size(v_detected,1), match.det_indices);
                    consistent = 0;

                    for extra = det_extra
                        v_extra = v_detected(extra,:);
                        v_rot = (R * v_extra')';

                        dist = pdist2(v_rot, v_catalog);
                        if min(dist) < 0.01
                            consistent = consistent + 1;
                        end
                    end

                    if consistent >= 2
                        matches{end+1} = match;
                        count = count + 1;

                        if count >= max_matches
                            return;
                        end
                    end
                end
            end
        end
    end
end

end