clc;
clear all;
close all;
addpath(genpath('cho_code'));
%%
weight_ring = 1;
%============================
%% Example 1
imgname = 'data/7_patch_use.png';
kernelname = [imgname(1:end-4) '_kernel.png'];
resultname = [imgname(1:end-4) '_out.png'];
%% parameter setting
lambda_tv = 0.01;
lambda_l0 = 2e-3;
%============================
% %% Example 2
% imgname = 'data\shelf.jpg';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.002;
% lambda_l0 = 2e-3;
% %============================
% %% Example 3
% imgname = 'data\8_patch_use.png';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.002;
% lambda_l0 = 2e-4;
% % %============================
% %% Example 4
% imgname = 'data\1_pattt_use.png';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.008;
% lambda_l0 = 2e-3;
% % %============================
% %% Example 5
% imgname = 'data\new_test_img12_blur_129.png';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.001;
% lambda_l0 = 2e-4;
% % %============================
% %% Example 6
% imgname = 'data\my_test_car6.png';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.002;
% lambda_l0 = 2e-3;
% % %============================
% %% Example 7
% imgname = 'data\eccv3_blurred.png';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.003;
% lambda_l0 = 1e-3;
% % %============================
% %% Example 8
% imgname = 'data\test3_blur_55.png';
% kernelname = [imgname(1:end-4) '_kernel.png'];
% resultname = [imgname(1:end-4) '_out.png'];
% %% parameter setting
% lambda_tv = 0.001;
% lambda_l0 = 5e-4;
% %============================
y = imread(imgname);
y = im2double(y);
k = imread(kernelname);
k = double(k(:,:,1));
kernel = k./sum(k(:));
%% Test examples
result = ringing_artifacts_removal(y, kernel, lambda_tv, lambda_l0, weight_ring);
figure; imshow(result,[]);
imwrite(result,resultname);

