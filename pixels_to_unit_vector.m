function [ux, uy, uz] = pixels_to_unit_vector(x_pix, y_pix)
% pixels_to_unit_vector - Convert undistorted pixel coordinates to unit vectors
%
%   Inputs:
%       xu_pix, yu_pix - undistorted pixel coordinates (scalars or vectors)
%
%   Outputs:
%       ux, uy, uz - components of the corresponding unit vector(s)
%
%   Uses Equation 4.1.6 (inverse pinhole model)

% --- camera parameters (same as distortion_correction.m) ---
f_mm = 25;                  % focal length [mm]
sensor_w = 5.70;            % sensor width [mm]
sensor_h = 4.28;            % sensor height [mm]
img_w = 646;                % image width [px]
img_h = 486;                % image height [px]

% pixel pitch (mm/px)
ppx = sensor_w / img_w;
ppy = sensor_h / img_h;

% principal point (assume at center)
xc = img_w / 2;
yc = img_h / 2;

% --- normalized coordinates (dimensionless) ---
alpha = (x_pix - xc) * ppx / f_mm; % x = column
beta  = (y_pix - yc) * ppy / f_mm; % y = row

% --- normalization factor ---
denom = sqrt(1 + alpha.^2 + beta.^2);

% --- unit vector components ---
ux = alpha ./ denom;
uy = beta  ./ denom;
uz = 1     ./ denom;
end