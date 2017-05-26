clc;
clear;
close all;
addpath(genpath('image'));
%parameter setting
opt.lambda_kernel = 0.01;% 0.01
opt.lambda_smooth = 0.005;%
opt.lambda_texture = 1;
opt.window_size = 45;
opt.lambda_kernel_smooth = 1e-5; %(Adjustable, typically 1e-5, 1e-4, 1e-3)
% if the kernel contains obvious noise, large value may help (e.g., for test1_blur.png (1e-3))
opt.display = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% example 1
% filename = 'image/details_blur2_45.png';
% filename = 'image/image2053R1_new.bmp';
% filename = 'image/wall.png';
% filename = 'test1_blur.png';
% filename = 'image/roma.png';
filename = 'image/image9.png';
% filename = 'image/grass_blur_using_Levin01_04_27.png';
% filename = 'image/structure.bmp';
% filename = 'lyndsey.png';
% filename = 'image10_large.png';
% filename ='chalk_drawing.jpg';

opt.kernel_sizeh = 35; opt.kernel_sizew = 35;

% kernel_sizeh = 35; kernel_sizew = 35;

opt.lambda_kernel_smooth = 1e-5;
% %% example 2
% filename = 'image\grass_blur_using_Levin01_04_27.png';
% kernel_sizeh = 27; kernel_sizew = 27;
% lambda_kernel_smooth = 1e-5;
% %% example 3
% filename = 'image\image2053R1_new.bmp';
% lambda_kernel_smooth = 1e-3;
% kernel_sizeh = 53; kernel_sizew = 53;
% %% example 4
% filename = 'image\lyndsey.png';
% lambda_kernel_smooth = 1e-4;
% kernel_sizeh = 31; kernel_sizew = 31;
% %% example 5
% filename = 'image\test1_blur.png';
% lambda_kernel_smooth = 1e-3;
% kernel_sizeh = 29; kernel_sizew = 29;
% %% example 6
% filename = 'image\roma.png';
% lambda_kernel_smooth = 1e-4;
% kernel_sizeh = 55; kernel_sizew = 55;
% %% example 7
% filename = 'image\wall.png';
% lambda_kernel_smooth = 1e-4;
% kernel_sizeh = 55; kernel_sizew = 55;
% %% example 8
% filename = 'image\cvpr_noise_6.bmp';
% lambda_kernel_smooth = 1e-4;
% kernel_sizeh = 35; kernel_sizew = 35;
% %% example 9
% filename = 'image\postcard.png';
% lambda_kernel_smooth = 1e-3;
% kernel_sizeh = 101; kernel_sizew = 57;
% %% example 10
% filename = 'image\structure.bmp';
% lambda_kernel_smooth = 1e-3;
% kernel_sizeh = 45; kernel_sizew = 45;
%============================

y = imread(filename);
if size(y,3)==3
    y_gray = rgb2gray(y);
    y_gray = double(y_gray)/255;
else
    y_gray = double(y)/255;
end

opt.kernel_size = 25;
times = 16;
patch_size = opt.kernel_size * times;

if patch_size < size(y_gray,1) && patch_size < size(y_gray,2)
    [pxmin,pxmax,pymin,pymax] = patchSelection(filename,opt.kernel_size,times);
    y_gray_P = y_gray(pxmin:pxmax,pymin:pymax);
else
    y_gray_P = y_gray;
end

% y_gray_P = imresize(y_gray,0.5);
ori_B =  im2double(y);

tic
[kernel, inter_latent] = deconv_main(y_gray_P, ori_B, opt);
toc

figure, imagesc(kernel)
figure, imshow(inter_latent)


% Imgoutname = [filename(1:end-4) '_out.png'];
% Kerneloutname = [filename(1:end-4) '_kernel.png'];
% kw= kernel-min(kernel(:));
% kw= kw/max(kw(:));
% imwrite(kw,Kerneloutname);
% imwrite(inter_latent,Imgoutname);

    
    
    
    
    
    
    
    