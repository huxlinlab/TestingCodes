function [fx, fy, rho, theta] = freq_coords(im_size)
%FREQ_COORDS get frequency domain (freq_x, freq_y) coordinates for 2d
%fourier domain. Takes care of tricky indexing issues for even/odd image
%sizes.
%
% [fx, fy] = FREQ_COORDS(im_size) returns meshgrid of coordinates fx and fy
% with the origin (DC component) in the center.
%
% [fx, fy, rho, theta] = FREQ_COORDS(im_size) includes meshgrid of rho
% (spatial frequency) and theta (angle), with theta=0 corresponding to
% vertical and increasing theta corresponding to CCW rotations. 'rho' has
% units of 'cycles', not 'cycles per pixel'.
%
% NOTE 1: theta is NaN at fx = 0, fy = 0
%
% NOTE 2: ifft2 expects the DC component at (1,1). Fourier filters derived
% from these [fx, fy, rho, theta] must be passed through ifftshift! For
% example:
%
%    noiseF = fft2(randn(im_size));
%    [~, ~, rho, ~] = FREQ_COORDS(im_size);
%    lowpass = exp(-rho.^2);
%    noiseF = noiseF .* lowpass;
%    img = real(ifft2(ifftshift(noiseF)));

% Note on even/odd image sizes: the DC component will be in the center for
% odd sizes, and at the ceiling of halfway for even sizes.
%
% Odd:  -4 -3 -2 -1 0 +1 +2 +3 +4
% Even: -4 -3 -2 -1 0 +1 +2 +3

range_start = -ceil((im_size - 1) / 2);
range_end = range_start + im_size - 1;

[fx, fy] = meshgrid(range_start:range_end);

if nargout >= 3
    rho = sqrt(fx.^2 + fy.^2);
end

if nargout >= 4
    theta = atan2(fx, -fy);
end
end