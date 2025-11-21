function [im, imF, filterF, aperture] = genImagesDT(width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa,aperture)
% Create a sequence band-pass grating (bpg) width x width images.
% spFreqCPP sets the mean % spatial frequency in cycles per pixel. 
% spFreqStdCPP sets the range of spatial frequencies present. 
% oriDEG sets the mean rotation
% oriKappa sets the range of orientation energy present.

noise = randn(width, width);
noiseF = fftshift(fft2(noise));%,width,width));
[~, ~, rho, theta] = freq_coords(width);

% Create separate [-1, 1] range meshgrid for pixel-space filters.
[px, py] = meshgrid(linspace(-1, 1, width));
pr = sqrt(px.^2 + py.^2);
im = zeros(width, width);imF = zeros(width, width);

% Create spatial frequency filter and gaussian aperture
spFreqFilter = pdf('rician', rho / width, spFreqCPP, spFreqStdCPP);
%aperture = exp(-4 * pr.^2);
%aperture = aperture*0+1;
% Create orientation filters. Note that 'theta' is doubled to create 
%  two symmetric filters in the Fourier domain (bow-tie rather
%  than cone shape). 'oriDEG' must also be doubled to compensate.
oriFilter = exp(oriKappa * cos((2 * theta) - (2 * deg2rad(oriDEG)))) / besseli(0, oriKappa);

% Get full, normalized foureir-domain filter.
filterF = spFreqFilter .* oriFilter;
filterF = filterF / sum(filterF(:));

% Apply fourier-domain filters on each frame.
imF(:, :) = squeeze(noiseF(:, :)) .* filterF;
im(:, :) = aperture .* real(ifft2(ifftshift(squeeze(imF(:, :)))));


% Normalize range in pixel space to +/- 1
im = im / max(abs(im(:)));
end

