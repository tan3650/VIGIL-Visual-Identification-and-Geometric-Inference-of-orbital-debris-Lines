img = imread('C:\Users\Test\OneDrive\Documents\Ultimate_Star\STARTRACKER\images\bmp_4.bmp');
img_double = im2double(img);     % Original image (for centroiding)
img_copy   = img_double;         % Copy for region growing (pixels get zeroed)

threshold  = 0.4; 
min_pixels = 4; 
max_pixels = 50; 

subset = img_copy(4:4:end, 2:2:end);
logical_subset = subset > threshold;

figure; imshow(img_double, []); hold on; % One figure for everything

if any(logical_subset(:))
    fprintf('At least one star candidate detected.\n');
    [rows, cols] = find(logical_subset);
    for k = 1:length(rows)
        r = rows(k); c = cols(k);
        full_r = r*4; full_c = c*2;
        fprintf('Candidate at (%d, %d)\n', full_r, full_c);

        [star_region, img_copy] = region_growing(img_copy, full_r, full_c, threshold, min_pixels, max_pixels, 0);

        if ~isempty(star_region)
            fprintf(' Valid star: %d px\n', size(star_region,1));
            [cx, cy] = centroiding(img_double, star_region, 5);
            fprintf(' Centroid: (%.2f, %.2f)\n', cy,cx);

            % Plot region + centroid
            plot(star_region(:,2), star_region(:,1), 'g.', 'MarkerSize', 6);
            plot(cy,cx, 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
        else
            fprintf(' Rejected.\n');
        end
    end
else
    fprintf('No candidates detected.\n');
end