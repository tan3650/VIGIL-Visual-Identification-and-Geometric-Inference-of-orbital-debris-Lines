function [cx_u, cy_u] = centroiding(img, star_region, window_size)
% CENTROIDING  Compute sub-pixel centroid of a detected star and correct for distortion
%
%   [cx_u, cy_u] = centroiding(img, star_region, window_size)
%
%   Inputs:
%       img         - 2D image matrix
%       star_region - Nx2 list of [row, col] pixels belonging to the star
%       window_size - size of square window around brightest pixel (default = 5)
%
%   Outputs:
%       cx_u - undistorted centroid row-position (vertical, sub-pixel)
%       cy_u - undistorted centroid column-position (horizontal, sub-pixel)

if nargin < 3
    window_size = 5; % default 5x5 window
end
half_win = floor(window_size / 2);

[rows, cols] = size(img);

% Find brightest pixel in region
intensities = zeros(size(star_region,1),1);
for k = 1:size(star_region,1)
    intensities(k) = img(star_region(k,1), star_region(k,2));
end
[~, idx_max] = max(intensities);
brightest_r = star_region(idx_max,1);
brightest_c = star_region(idx_max,2);

% Define window around brightest pixel
r_min = max(1, brightest_r - half_win);
r_max = min(rows, brightest_r + half_win);
c_min = max(1, brightest_c - half_win);
c_max = min(cols, brightest_c + half_win);

[win_cols, win_rows] = meshgrid(c_min:c_max, r_min:r_max);
win_rows = win_rows(:);
win_cols = win_cols(:);

% Grab intensities from the window
win_intensities = zeros(length(win_rows),1);
for k = 1:length(win_rows)
    win_intensities(k) = img(win_rows(k), win_cols(k));
end

% Weighted centroid (center of gravity)
total_intensity = sum(win_intensities);
if total_intensity == 0
    cx_d = NaN;
    cy_d = NaN;
else
    cy_d = sum(win_cols .* win_intensities) / total_intensity; % distorted col
    cx_d = sum(win_rows .* win_intensities) / total_intensity; % distorted row
end

% ---- Apply distortion correction ----
if ~isnan(cx_d) && ~isnan(cy_d)

    % cx_d = row (y), cy_d = col (x)
    x_d = cy_d;   % column
    y_d = cx_d;   % row

    % distortion_correction expects (x, y)
    [x_u, y_u] = distortion_correction(x_d, y_d);

    % return in SAME format as before (row, col)
    cx_u = y_u;   % row
    cy_u = x_u;   % col

    fprintf('Centroid (distorted):   row = %.3f, col = %.3f\n', cx_d, cy_d);
    fprintf('Centroid (undistorted): row = %.3f, col = %.3f\n', cx_u, cy_u);

else
    cx_u = NaN;
    cy_u = NaN;
    fprintf('Centroid not found (zero intensity).\n');
end
end