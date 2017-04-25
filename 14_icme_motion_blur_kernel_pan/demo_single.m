clc;
clear;
close all;
addpath(genpath('image'));
addpath(genpath('denoising'));
addpath(genpath('Levin_deblurring_code'));
 %parameter setting
 lambda_kernel = 0.01;% 0.01
 lambda_smooth = 0.005;%
 % lambda_texture = 0.02; %for L0Smoothing
 lambda_texture = 1;
 window_size = 3;
 gamma_correction = 1.0;
 lambda_kernel_smooth = 1e-4;
 
display = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
y = imread('image/im01_ker01_blur.png');
% Very important paramters
kernel_sizeh =19;
kernel_sizew =19;
%
% %%
% y_im = structure_adaptive_map(im2double(y), 1/15, 100, 0);
% y_im = y_im-min(y_im(:));
% y_im = y_im/max(y_im(:));
% imwrite(y_im, 'noise_image\dlut_blur_denoise.png')
% %%
% y = imread('noise_image\dlut_blur_denoise.png');
% % %%%%%%%%%%%%
if (size(y,3) ==3)
    B_patch = rgb2gray(y);
    B_patch = double(B_patch)/255;
else
    B_patch = double(y)/255;
end
ori_B =  im2double(y);
B_patch = B_patch.^gamma_correction;
%[kernel_sizeh, kernel_sizew]= size(f);
tic;
[deblur, kernel] = deconv_main(B_patch, ori_B, lambda_kernel, lambda_smooth,...
    lambda_texture, window_size, kernel_sizeh, kernel_sizew, ...
    lambda_kernel_smooth, display);
toc
% for i = 1:size(ori_B,3)
%     deblur_f(:,:,i) = deblurring_lr_whole(ori_B(:,:,i), kernel, 0.01);
% end
%% Final Deblurring
% lambda = 4e3;
% alpha = 2/3;
% %% Boundary processing
% H = size(ori_B,1);    W = size(ori_B,2);
% y_pad = wrap_boundary_liu(ori_B, opt_fft_size([H W]+size(kernel)-1));
% for c = size(ori_B,3)
%     deblur_f(:,:,c) = fast_deconv(y_pad(:,:,c), flp(kernel), lambda, alpha);
% end
% deblur_f = deblur_f(1:H,1:W,:);
% figure; imshow(deblur_f)
lambda = 0.003;
clear deblur_f;
for i = 1:size(ori_B,3)
    deblur_f(:,:,i) = deconvSps(ori_B(:,:,i), kernel, lambda, 200);
end
%save('Blurry2_4','kernel','deblur_f');
    
    
    
    
    
    
    
    