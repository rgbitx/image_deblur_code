clc;
clear all;
close all;
addpath(genpath('image'));
addpath(genpath('whyte_code'));
addpath(genpath('cho_code'));
opts.prescale = 1; %%downsampling
opts.xk_iter = 5; %% the iterations
opts.gamma_correct = 1.0;
opts.k_thresh = 20;
%% Test on 2013-08-17
% filename = 'image\test3_blur_55.png'; opts.kernel_size = 69;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% %% 
% filename = 'image\0015_blur65.png'; opts.kernel_size =67;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-4; weight_ring = 1;
%
% filename = 'image\7_patch_use.png'; opts.kernel_size = 85;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; 
% lambda_tv = 0.01; lambda_l0 = 2e-3; weight_ring = 1;
%
% filename = 'image\1_pattt_use.png'; opts.kernel_size = 55;   saturation = 0;
% lambda_pixel = 0; lambda_grad = 4e-3; opts.gamma_correct = 2.2;
% lambda_tv = 0.008; lambda_l0 = 2e-3; weight_ring = 1;
%
filename = 'image/8_patch_use.png'; opts.kernel_size = 25;  saturation = 0;
lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 2.2;
lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
%
% filename = 'image\dragon_patch_use.png'; opts.kernel_size = 45;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.003; lambda_l0 = 2e-3; weight_ring = 1;
%
% filename = 'image\new_test_img_blur.png'; opts.kernel_size = 55;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 1e-4; weight_ring = 1;
% %%
% filename = 'image\2013622235456945_blur_79.png'; opts.kernel_size = 85;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0001; lambda_l0 = 1e-4; weight_ring = 1;
%%
% filename = 'image\my_shufa_image_blur_121.png'; opts.kernel_size = 111;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0002; lambda_l0 = 1e-4; weight_ring = 1;
%
% filename = 'image\new_test_img10_blur.png'; opts.kernel_size = 85; saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0002; lambda_l0 = 1e-4;
%
% filename = 'image\new_test_img12_blur_129.png'; opts.kernel_size = 131;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-4; weight_ring = 1;
%
% filename = 'image\eccv3_blurred.png'; opts.kernel_size = 23;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.003; lambda_l0 = 1e-3; weight_ring = 1;
%%
% filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
%%
% filename = 'image\color_patch_blur_99.png'; opts.kernel_size = 99;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% %%
% filename = 'image\blurred_cho.png'; opts.kernel_size = 23;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.0002; lambda_l0 = 2e-4; weight_ring = 1;
%%
% filename = 'image\my_test_car6.png'; opts.kernel_size = 95;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
%%
% filename = 'image\IMG_0747_Gaussion.png'; opts.kernel_size = 75;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
%%
% filename = 'image\DSC0065_small.png'; opts.kernel_size = 99;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% %%
% filename = 'image\image_saturation_1.png'; opts.kernel_size = 69;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
%
% filename = 'image\blurred.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
%%
% filename = 'image\IMG_0657_sss.JPG'; opts.kernel_size = 35; saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% %%
% filename = 'image\IMG_0141_small.jpg'; opts.kernel_size = 53; saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
%%
% filename = 'image\street_cars_blurry.jpg'; opts.kernel_size = 53;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
%%
%% Natual images deblurring examples
% filename = 'image\boat_input.png'; opts.kernel_size = 25;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
%%
% filename = 'image\hanzi.jpg'; opts.kernel_size = 33;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.003; lambda_l0 = 2e-3; weight_ring = 1;
%%
% filename = 'image\fountain1_blurry.png'; opts.kernel_size = 41;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 0;
%%
% filename = 'image\postcard.png'; opts.kernel_size = 91;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0005; lambda_l0 = 2e-3; weight_ring = 0;
%=================================
% filename = 'image\im06_ker01_blur.png'; opts.kernel_size = 19;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 0;
% filename = 'image\im15_ker04_blur.png'; opts.kernel_size = 27;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 0;
%===================================
y = imread(filename);
% y = y(3:end-2,3:end-2,:);
%y = imfilter(y,fspecial('gaussian',5,2),'same','replicate'); 
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
[kernel, interim_latent] = blind_deconv(yg, lambda_pixel, lambda_grad, opts);
toc
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

