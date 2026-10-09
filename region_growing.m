function [star_region, img] = region_growing(img, seed_r, seed_c, threshold, min_pixels, max_pixels, show_flag)
% REGION_GROWING  Iterative region growing for star detection + optional visualization
%
%   [star_region, img] = region_growing(img, seed_r, seed_c, threshold, min_pixels, max_pixels, show_flag)
%
%   Inputs:
%       img        - 2D image matrix (double or uint16 etc.)
%       seed_r     - seed row index
%       seed_c     - seed column index
%       threshold  - pixel intensity threshold
%       min_pixels - minimum allowed pixels in a star
%       max_pixels - maximum allowed pixels in a star
%       show_flag  - set to true (1) to display region overlay
%
%   Output:
%       star_region - Nx2 matrix of [row, col] pixel indices belonging to the star
%       img         - modified image with detected star pixels zeroed

[rows, cols] = size(img);
stack = [seed_r, seed_c];
region_pixels = [];

while ~isempty(stack)
    pixel = stack(end, :);
    stack(end, :) = [];
    r = pixel(1); c = pixel(2);
    if r < 1 || r > rows || c < 1 || c > cols
        continue;
    end
    if img(r, c) > threshold
        region_pixels(end+1, :) = [r, c]; %#ok<AGROW>
        img(r, c) = 0;
        stack(end+1, :) = [r-1, c]; %#ok<AGROW>
        stack(end+1, :) = [r+1, c]; %#ok<AGROW>
        stack(end+1, :) = [r, c-1]; %#ok<AGROW>
        stack(end+1, :) = [r, c+1]; %#ok<AGROW>
    end
end

n_pixels = size(region_pixels, 1);
if n_pixels < min_pixels || n_pixels > max_pixels
    star_region = [];
else
    star_region = region_pixels;
end

end