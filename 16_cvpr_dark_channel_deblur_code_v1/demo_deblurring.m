clc;
clear;
close all;
addpath(genpath('image'));
addpath(genpath('whyte_code'));
addpath(genpath('cho_code'));
opts.prescale = 1; %%downsampling
opts.xk_iter = 5; %% the iterations
opts.gamma_correct = 1.0;
opts.k_thresh = 20;
%% Note:
%% lambda_tv, lambda_l0, weight_ring are non-necessary, they are not used in kernel estimation.
%%
% filename = 'image\7_patch_use.png'; opts.kernel_size = 85;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;
% lambda_tv = 0.01; lambda_l0 = 2e-3; weight_ring = 1;
%%
% filename = 'image\flower.jpg'; opts.kernel_size = 35;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3; 
% lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;
%%
% filename = 'image\summerhouse.jpg'; opts.kernel_size = 95;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3; 
% lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;
%%
% filename = 'image\postcard.png'; opts.kernel_size = 115;  saturation =0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\boat.jpg'; opts.kernel_size = 35;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\flower_blurred.png'; opts.kernel_size = 55;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\wall.png'; opts.kernel_size = 65;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0001; lambda_l0 = 2e-3; weight_ring = 0;
%%
% filename = 'image\blurry_2_small.png'; opts.kernel_size = 35;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\blurry_7.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_dark = 0; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\my_test_car6.png'; opts.kernel_size = 95;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% filename = 'BlurryImages\Blurry4_9.png'; opts.kernel_size = 99;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 1e-3; weight_ring = 0;
%
% filename = 'BlurryImages\Blurry4_6.png'; opts.kernel_size = 41;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 1e-3; weight_ring = 0;
% filename = 'image\toy.png'; opts.kernel_size = 101;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\Blurry2_10.png'; opts.kernel_size = 105;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\im05_ker04_blur.png'; opts.kernel_size = 27;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_1240_blur.png'; opts.kernel_size = 45;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\IMG_0650_small_patch.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\IMG_0664_small_patch.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%% For the first figure

% filename = 'image\real_leaffiltered.png'; opts.kernel_size = 65;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1; % Gaussian filters

%%
% filename = 'image\IMG_4548_small.png'; opts.kernel_size = 35;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1; % Gaussian filters
%%
% filename = 'image\las_vegas_saturated.png'; opts.kernel_size = 99;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\IMG_4561.JPG'; opts.kernel_size = 65;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\IMG_4355_small.png'; opts.kernel_size = 45;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
%
% filename = 'image\IMG_4528_patch.png'; opts.kernel_size = 75;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\IMG_4528_patch.png'; opts.kernel_size = 75;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\26.blurred.jpg'; opts.kernel_size = 45;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% %%
% filename = 'image\26.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
%%
% filename = 'image\real_blur_img3.png'; opts.kernel_size = 35;  saturation = 0;
% lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
filename = 'image\real_img2.png'; opts.kernel_size = 25;  saturation = 0;
lambda_dark = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
lambda_tv = 0.003; lambda_l0 = 5e-4; weight_ring = 1;
%%
%===================================
y = imread(filename);
% y = y(3:end-2,3:end-2,:);
% y = imfilter(y,fspecial('gaussian',5,1),'same','replicate'); 
isselect = 0; %false or true
if isselect ==1
    figure, imshow(y);
    %tips = msgbox('Please choose the area for deblurring:');
    fprintf('Please choose the area for deblurring:\n');
    h = imrect;
    position = wait(h);
    close;
    B_patch = imcrop(y,position);
    y = (B_patch);
else
    y = y;
end
if size(y,3)==3
    yg = im2double(rgb2gray(y));
else
    yg = im2double(y);
end
tic;
[kernel, interim_latent] = blind_deconv(yg, lambda_dark, lambda_grad, opts);
toc
%% Algorithm is done!
%% ============Non-blind deconvolution ((Just use text deblurring methods))====================%%
y = im2double(y);
%% Final Deblur: 
if ~saturation
    %% 1. TV-L2 denoising method
    Latent = ringing_artifacts_removal(y, kernel, lambda_tv, lambda_l0, weight_ring);
else
    %% 2. Whyte's deconvolution method (For saturated images)
    Latent = whyte_deconv(y, kernel);
end
figure; imshow(Latent)
%%
k = kernel - min(kernel(:));
k = k./max(k(:));
imwrite(k,['results\' filename(7:end-4) '_kernel.png']);
imwrite(Latent,['results\' filename(7:end-4) '_result.png']);
imwrite(interim_latent,['results\' filename(7:end-4) '_interim_result.png']);
%%
%% ============Use Cho et al., ICCV 11 for non-blind deconvolution (For saturated images???)====================%%
%deblurred0 = deconv_outlier(y, kernel, 5/255, 0.003);
