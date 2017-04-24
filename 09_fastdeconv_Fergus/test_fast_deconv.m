% clear
% Test the fast deconvolution  method presented in the paper
% D. Krishnan, R. Fergus: "Fast Image Deconvolution using Hyper-Laplacian
% Priors", Proceedings of NIPS 2009.

% parameter values; other values such as the continuation regime of the
% parameter beta should be changed in fast_deconv.m
lambda = 2e3;
alpha = 1/2;

% read in and normalize input data - if data is not normalized, change
% lambda parameter accordingly
y = double(imread('dsc_0085.jpg'))./255.0;

% load 2 pre-defined kernels
load kernels.mat;
% which kernel to use
kernel = kernel1;
ks = floor((size(kernel, 1) - 1)/2);

[imy,imx]=size(y);
hy = floor(imy/2); hx = floor(imx/2);
% crop to some pre-defined size
cy = 1872; cx = 1866;
% cy=1600;cx=5000;
IMG_SIZE = 1024;
% IMG_SIZE = 1024;
yorig = y(cy-IMG_SIZE/2-ks:cy+IMG_SIZE/2+ks, cx-IMG_SIZE/2-ks:cx+IMG_SIZE/2+ks);
y = y(cy-IMG_SIZE/2:cy+IMG_SIZE/2, cx-IMG_SIZE/2:cx+IMG_SIZE/2);

% convolve with kernel1 and add noise
yorig = y;
y = conv2(yorig, kernel1, 'valid');
y = y + 0.01*randn(size(y));
y = double(uint8(y .* 255))./255;

% edgetaper to better handle circular boundary conditions
y = padarray(y, [1 1]*ks, 'replicate', 'both');
% for a=1:4
%   y = edgetaper(y, kernel);
% end

% Check if Eero Simoncell's function exists
if (exist('pointOp') ~= 3) 
  fprintf('WARNING: Will use slower interp1 for LUT interpolations. For speed please see comments at the top of fast_deconv.m\n'); 
end;

clear persistent;
snr_blur = snr1(y, ks, yorig);
tic;
[x] = fast_deconv(y, kernel, lambda, alpha);
timetaken = toc;
snr_recon = snr1(x, ks);

% remove padding
x = x(ks+1:end-ks, ks+1:end-ks);
y = y(ks+1:end-ks, ks+1:end-ks);
yorig = yorig(ks+1:end-ks, ks+1:end-ks);

figure; imagesc([yorig y x]); colormap gray; 
tt = sprintf('Original  Blurred (SNR %.2f) Reconstructed (SNR %.2f)', ...
             snr_blur, snr_recon);
title(tt);

fprintf('Time taken for image of size %dx%d is %.3f\n', IMG_SIZE, IMG_SIZE, ...
        timetaken);

