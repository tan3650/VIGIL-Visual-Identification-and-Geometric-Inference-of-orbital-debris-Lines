clear; clc; close all;

%% --- Step 1: Load Image + inject debris and noise ---
orig = im2double(imread('images/bmp_4.bmp'));
if size(orig,3)==3, orig = rgb2gray(orig); end
[H,W] = size(orig);
[X,Y] = meshgrid(1:W, 1:H);

% adaptive star mask (stars are restored untouched at the end)
bg    = median(orig(:));
sigma = 1.4826*median(abs(orig(:)-bg));
thr   = max(bg + 8*max(sigma,0.005), 0.15);
raw   = bwareaopen(orig > thr, 4);
starMask = imdilate(raw, strel('disk',1));
fprintf('bg=%.3f thr=%.3f star mask covers %.2f%% of pixels\n', bg, thr, 100*mean(starMask(:)));

% debris streaks (own layer)
debris = zeros(H,W);
nStreaks = 6;
truth = zeros(nStreaks,5);                             % x0 y0 x1 y1 amp
baseAng = 2*pi*rand;                                   % shared direction
for k = 1:nStreaks
    L   = round(W*(0.10 + 0.15*rand));
    ang = baseAng + 0.25*randn;                        % similar directions
    x0  = W*(0.1+0.8*rand);  y0 = H*(0.1+0.8*rand);
    x1  = x0 + L*cos(ang);   y1 = y0 + L*sin(ang);
    amp = 0.06 + 0.14*rand;                            % faint
    sig = 1.0 + 0.8*rand;
    dx = x1-x0;  dy = y1-y0;
    t  = ((X-x0)*dx + (Y-y0)*dy)/(dx^2+dy^2);
    tc = min(max(t,0),1);
    d2 = (X-(x0+tc*dx)).^2 + (Y-(y0+tc*dy)).^2;
    taper = min(1, 6*min(tc,1-tc)) .* (t>=-0.02 & t<=1.02);   % soft ends
    flash = 1 + 0.3*sin(2*pi*(3*tc + rand));                  % tumble
    debris = debris + amp*flash.*taper.*exp(-d2/(2*sig^2));
    truth(k,:) = [x0 y0 x1 y1 amp];
end

% build on the ORIGINAL image: debris + sky + noise
layer = orig + debris + 0.03 + 0.03*(X/W);
layer = layer + sqrt(max(layer,0)/1000).*randn(H,W);   % shot noise
layer = layer + 0.01*randn(H,W);                       % read noise

s = rng; rng(1);                                       % fixed hot/dead pixels
layer(randperm(H*W, round(0.0005*H*W))) = 1;
layer(randperm(H*W, round(0.0003*H*W))) = 0;
rng(s);

for k = 1:15                                           % cosmic rays
    cx = randi(W); cy = randi(H);
    for j = 0:randi(3)-1
        layer(cy, min(cx+j,W)) = 0.6 + 0.4*rand;
    end
end
layer = imnoise(min(max(layer,0),1), 'salt & pepper', 0.002);

img = layer;
img(starMask) = orig(starMask);                        % stars untouched
imwrite(img, 'images/bmp_4_noisy.png');

figure; imshow(img); title('Step 1: Raw Input Image');

%% --- Step 2: Detection (stars vs debris) ---
threshold = 0.4; min_pixels = 4; max_pixels = 50; window_size = 5;
img_double = im2double(img);

bw = bwareaopen(img_double > threshold, min_pixels);
cc = bwconncomp(bw);
stats = regionprops(cc, 'Area','MajorAxisLength','MinorAxisLength');
elong  = [stats.MajorAxisLength] ./ max([stats.MinorAxisLength], 1);
isStar = ([stats.Area] <= max_pixels) & (elong < 3);

lbl    = labelmatrix(cc);
starBW = ismember(lbl, find(isStar));
fprintf('Star components: %d | Large/elongated components: %d (injected streaks: %d)\n', ...
        sum(isStar), sum(~isStar), nStreaks);

img_det  = img_double .* starBW;     % only star-like objects survive
img_copy = img_det;

subset = img_copy(4:4:end, 2:2:end);
logical_subset = subset > threshold;

v_detected = [];
pixels = [];

if any(logical_subset(:))
    [rows, cols] = find(logical_subset);
    figure; imshow(img_double, []); hold on;

    for k = 1:length(rows)
        r = rows(k); c = cols(k);
        full_r = r*4; full_c = c*2;

        [star_region, img_copy] = region_growing(img_copy, full_r, full_c, threshold, min_pixels, max_pixels, 0);

        if ~isempty(star_region)
            [cx_u, cy_u] = centroiding(img_det, star_region, window_size);

            x = cy_u; y = cx_u;
            [ux, uy, uz] = pixels_to_unit_vector(x, y);

            v_detected = [v_detected; ux, uy, uz];
            pixels = [pixels; cy_u, cx_u];

            plot(star_region(:,2), star_region(:,1), 'g.');
            plot(cy_u, cx_u, 'rx');
        end
    end
else
    error('No stars detected.');
end

title('Step 2: Detected Stars');
hold off;

%% --- Step A: Detected vectors ---
figure;
[X,Y,Z] = sphere(50);
surf(X,Y,Z,'FaceAlpha',0.05,'EdgeColor','none'); hold on;
axis equal; grid on;

for i = 1:size(v_detected,1)
    quiver3(0,0,0, v_detected(i,1), v_detected(i,2), v_detected(i,3),'r','LineWidth',2);
end

title('Step A: Detected Vectors (Camera Frame)');

%% --- Step 3: Catalog ---
[v_catalog, star_names, ra_deg, dec_deg, mag_catalog] = ...
    catalog_matching('catalogues/gemini.csv');

%% --- Step B: Catalog vectors ---
figure;
[X,Y,Z] = sphere(50);
surf(X,Y,Z,'FaceAlpha',0.05,'EdgeColor','none'); hold on;
axis equal; grid on;

plot3(v_catalog(:,1), v_catalog(:,2), v_catalog(:,3), 'b.');

title('Step B: Catalog Vectors');

%% --- Step 4: DB ---
dbtol = 0.1;
triangle_db = build_triangle_db(v_catalog, dbtol);

%% --- Step 5: Matching ---
matches = pattern_matching(v_detected, triangle_db, v_catalog, dbtol);
if isempty(matches), error('No matches found'); end

%% --- Step 6: Attitude ---
[R_best, bestMatch, v_rotated] = ...
    attitude_determination(v_detected, v_catalog, matches);

%% --- Step C: Triangle match ---
figure;
[X,Y,Z] = sphere(50);
surf(X,Y,Z,'FaceAlpha',0.05,'EdgeColor','none'); hold on;
axis equal; grid on;

v_det_tri = v_detected(bestMatch.det_indices,:);
v_cat_tri = v_catalog(bestMatch.cat_indices,:);

for i = 1:3
    quiver3(0,0,0, v_det_tri(i,1), v_det_tri(i,2), v_det_tri(i,3),'r','LineWidth',3);
    quiver3(0,0,0, v_cat_tri(i,1), v_cat_tri(i,2), v_cat_tri(i,3),'g','LineWidth',3);
end

title('Step C: Triangle Match (Red vs Green)');

%% --- Step D: Rotation ---
figure;
[X,Y,Z] = sphere(50);
surf(X,Y,Z,'FaceAlpha',0.05,'EdgeColor','none'); hold on;
axis equal; grid on;

for i = 1:size(v_detected,1)
    quiver3(0,0,0, v_detected(i,1), v_detected(i,2), v_detected(i,3),'r');
end

for i = 1:size(v_rotated,1)
    quiver3(0,0,0, v_rotated(i,1), v_rotated(i,2), v_rotated(i,3),'b','LineWidth',2);
end

title('Step D: Rotation (Red = Before, Blue = After)');

%% --- FINAL GLOBAL ASSIGNMENT ---
dist = pdist2(v_rotated, v_catalog);

final_indices = zeros(size(v_detected,1),1);
used = false(size(v_catalog,1),1);

for i = 1:size(v_detected,1)
    [~, idx] = sort(dist(i,:));

    for j = 1:length(idx)
        if ~used(idx(j))
            final_indices(i) = idx(j);
            used(idx(j)) = true;
            break;
        end
    end
end

%% -------- FINAL TABLE --------
N = size(v_detected,1);

ResultTable = table(...
    (1:N)', ...
    pixels(:,1), pixels(:,2), ...
    v_detected(:,1), v_detected(:,2), v_detected(:,3), ...
    v_rotated(:,1), v_rotated(:,2), v_rotated(:,3), ...
    v_catalog(final_indices,1), ...
    v_catalog(final_indices,2), ...
    v_catalog(final_indices,3), ...
    string(star_names(final_indices)), ...
    'VariableNames', { ...
    'StarID','Pixel_X','Pixel_Y', ...
    'Det_X','Det_Y','Det_Z', ...
    'Rot_X','Rot_Y','Rot_Z', ...
    'Cat_X','Cat_Y','Cat_Z', ...
    'MatchedStar'});

disp('================ FINAL STAR TABLE ================');
disp(ResultTable);

%% -------- FINAL ATTITUDE --------
disp('================ FINAL ATTITUDE ================');
disp('Rotation Matrix:');
disp(R_best);

yaw   = atan2d(R_best(2,1), R_best(1,1));
pitch = -asind(R_best(3,1));
roll  = atan2d(R_best(3,2), R_best(3,3));

fprintf('Euler Angles (deg):\n');
fprintf('Yaw   = %.3f\n', yaw);
fprintf('Pitch = %.3f\n', pitch);
fprintf('Roll  = %.3f\n', roll);

%% --- Step E: Final matching ---
figure;
[X,Y,Z] = sphere(50);
surf(X,Y,Z,'FaceAlpha',0.05,'EdgeColor','none'); hold on;
axis equal; grid on;

for i = 1:length(final_indices)
    cat_idx = final_indices(i);

    quiver3(0,0,0, v_rotated(i,1), v_rotated(i,2), v_rotated(i,3),'b','LineWidth',2);
    quiver3(0,0,0, v_catalog(cat_idx,1), v_catalog(cat_idx,2), v_catalog(cat_idx,3),'g','LineWidth',2);

    plot3([v_rotated(i,1) v_catalog(cat_idx,1)], ...
          [v_rotated(i,2) v_catalog(cat_idx,2)], ...
          [v_rotated(i,3) v_catalog(cat_idx,3)], 'k--');

    text(v_catalog(cat_idx,1), v_catalog(cat_idx,2), v_catalog(cat_idx,3), ...
        star_names{cat_idx}, 'Color','yellow');
end

title('Step E: Final Matching');

%% --- Image overlay ---
figure;
imshow(img, []); hold on;

plot(pixels(:,1), pixels(:,2), 'go','LineWidth',1.5);

for i = 1:length(final_indices)
    idx = final_indices(i);
    text(pixels(i,1)+10, pixels(i,2), star_names{idx}, ...
        'Color','yellow','FontWeight','bold');
end

title('Final Result on Image');
hold off;