function [xu_pix, yu_pix] = distortion_correction(xd_pix, yd_pix)
% distortion_correction - Corrects lens distortion using Brown's model
%
%   Inputs:
%       xd_pix, yd_pix - distorted pixel coordinates (can be scalars or vectors)
%
%   Outputs:
%       xu_pix, yu_pix - undistorted pixel coordinates
%
%   Based on fitted coefficients:
%       K1 = 1.1999879539599927e-07
%       K2 = -2.5811362783036343e-13
%       P1 = 6.604606889315861e-06
%       P2 = 1.5756800510123402e-05
%
%   Camera geometry:
%       f = 25 mm
%       sensor = 5.70 mm x 4.28 mm
%       image = 646 x 486 pixels
%

% coefficients
K1 = 1.1999879539599927e-07;
K2 = -2.5811362783036343e-13;
P1 = 6.604606889315861e-06;
P2 = 1.5756800510123402e-05;

% camera parameters
f = 25;                     % focal length [mm]
sensor_w = 5.70;            % sensor width [mm]
sensor_h = 4.28;            % sensor height [mm]
img_w = 646;                % image width [px]
img_h = 486;                % image height [px]

% focal length in pixels
fx = f * (img_w / sensor_w);
fy = f * (img_h / sensor_h);
fpix = mean([fx, fy]);      % assume square pixels

% principal point (center of image)
xc = img_w / 2;
yc = img_h / 2;

% --- normalize distorted coordinates ---
x = (xd_pix - xc) / fpix;
y = (yd_pix - yc) / fpix;

% radius squared
r2 = x.^2 + y.^2;

% --- Brown distortion model ---
xu = x .* (1 + K1*r2 + K2*r2.^2) + P2*(r2 + 2*x.^2) + 2*P1.*x.*y;
yu = y .* (1 + K1*r2 + K2*r2.^2) + P1*(r2 + 2*y.^2) + 2*P2.*x.*y;

% --- back to pixel coordinates ---
xu_pix = fpix * xu + xc;
yu_pix = fpix * yu + yc;

end

