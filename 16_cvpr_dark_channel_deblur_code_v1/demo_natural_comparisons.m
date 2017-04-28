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
 %%
% filename = 'image\flower.jpg'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; 
% lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;

% filename = 'image\result_n_comparison.bmp'; opts.kernel_size = 115;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; 
% lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;
% filename = 'image\summerhouse.jpg'; opts.kernel_size = 95;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; 
% lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;
% filename = 'image\26.blurred.jpg'; opts.kernel_size = 41;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;
%%
% filename = 'image\Blurry2_8.png'; opts.kernel_size = 113;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 0;
%%
% filename = 'image\postcard.png'; opts.kernel_size = 115;  saturation =0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\boat.jpg'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\flower_blurred.png'; opts.kernel_size = 55;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\wall.png'; opts.kernel_size = 65;  saturation = 0;
% lambda_pixel = 0; lambda_grad = 8e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0001; lambda_l0 = 2e-3; weight_ring = 0;
% filename = 'image\05-Seungyong-Lee.bmp'; opts.kernel_size = 55;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
%%
% filename = 'image\DSC0065_small.png'; opts.kernel_size = 99;  saturation = 1;
% lambda_pixel = 0; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\blurry_2_small.png'; opts.kernel_size = 35;  saturation = 1;
% lambda_pixel = 0; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\blurry_7.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_pixel = 0; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\my_test_car6.png'; opts.kernel_size = 95;  saturation = 1;
% lambda_pixel = 0; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% filename = 'BlurryImages\Blurry4_9.png'; opts.kernel_size = 99;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 1e-3; weight_ring = 0;
%
% filename = 'BlurryImages\Blurry4_6.png'; opts.kernel_size = 41;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 1e-3; weight_ring = 0;
% filename = 'image\toy.png'; opts.kernel_size = 101;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\Blurry2_10.png'; opts.kernel_size = 105;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\im05_ker04_blur.png'; opts.kernel_size = 27;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_1240_blur.png'; opts.kernel_size = 45;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\Untitled.png'; opts.kernel_size = 25;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_1245.png'; opts.kernel_size = 39;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\-5.5um.tif'; opts.kernel_size = 11;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\blur_sat1.png'; opts.kernel_size = 117;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;

% filename = 'image\real_leaf.png'; opts.kernel_size = 65;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_4548_smallfiltered.png'; opts.kernel_size = 35;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\moto_car.png'; opts.kernel_size = 51;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\26.png'; opts.kernel_size = 65;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\flower_blurred.png'; opts.kernel_size = 55;  saturation = 0;
% lambda_pixel = 0; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\roma.jpg'; opts.kernel_size = 85;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\boat_input.png'; opts.kernel_size = 25;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\lyndsey.png'; opts.kernel_size = 27;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
%%
% filename = 'image\boat.jpg'; opts.kernel_size = 31;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\boy_statue.jpg'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\fishes.jpg'; opts.kernel_size = 45;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\hanzi.jpg'; opts.kernel_size = 33;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\building.jpg'; opts.kernel_size = 33;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.003; lambda_l0 = 5e-4; weight_ring = 0;
% filename = 'image\harubang.jpg'; opts.kernel_size = 61;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\picassoBlurImage.png'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0001; lambda_l0 = 5e-4; weight_ring = 0;
%%
% filename = 'image\mukta.jpg'; opts.kernel_size = 25;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 0;
% filename = 'image\summerhouse.jpg'; opts.kernel_size = 115;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 0;
% filename = 'image\bird.png'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 0;
% filename = 'image\nvBlurImage.png'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\stuatus.jpg'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\8_patch_use.png'; opts.kernel_size = 125;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\fountain1_blurry.png'; opts.kernel_size = 45;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_4530_patch.png'; opts.kernel_size = 75;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_4528_patch.png'; opts.kernel_size = 75;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\IMG_4532_patch.JPG'; opts.kernel_size = 65;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\car_from_mit.png'; opts.kernel_size = 95;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\car_from_mit2.png'; opts.kernel_size = 95;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 1.0;
% lambda_tv = 0.0005; lambda_l0 = 5e-4; weight_ring = 1;
% filename = 'image\my_shufa_image_blur_121.png'; opts.kernel_size = 111;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0002; lambda_l0 = 1e-4; weight_ring = 1;
% filename = 'image\0015_blur65.png'; opts.kernel_size =67;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.001; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\2013622235456945_blur_79.png'; opts.kernel_size = 85;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% lambda_tv = 0.0001; lambda_l0 = 1e-4; weight_ring = 1;
% filename = 'image\1_pattt_use.png'; opts.kernel_size = 55;   saturation = 0;
% lambda_pixel = 0; lambda_grad = 4e-3; opts.gamma_correct = 2.2;
% lambda_tv = 0.008; lambda_l0 = 2e-3; weight_ring = 1;
% filename = 'image\color_patch_blur_99.png'; opts.kernel_size = 99;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\my_test_car6.png'; opts.kernel_size = 95;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;opts.gamma_correct = 2.2;
% filename = 'image\blurry_7.png'; opts.kernel_size = 55;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\blurry_2_small.png'; opts.kernel_size = 35;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\IMG_3164_small.png'; opts.kernel_size = 45;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\real_img2.png'; opts.kernel_size = 29;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\030_03_01_051_05_large_blur.bmp'; opts.kernel_size = 75;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\face02_blurred.bmp'; opts.kernel_size = 25;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% filename = 'image\26.blurred.jpg'; opts.kernel_size = 45;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 0;
% filename = 'image\IMG_4350_small.png'; opts.kernel_size = 75;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 0;
filename = 'image\IMG_4355_small.png'; opts.kernel_size = 45;  saturation = 0;
lambda_pixel = 4e-3; lambda_grad = 4e-3;%opts.gamma_correct = 2.2;
lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
%===================================
y = imread(filename);
% y = y(3:end-2,3:end-2,:);
%y = imfilter(y,fspecial('gaussian',5,1),'same','replicate'); 
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
%deblurred0 = deconv_outlier(y, kernel, 12/255, 0.003);
