clear,clc,close all hidden;
% addpath('cho_code');
%% get input
%img_name_set = {'stone1_crop.png', 'stone2_crop.png'};
img_name_set = {'book1.png', 'book3.png'};%,'book3.png'};
% img_name_set = {'stone1.png', 'stone2.png'};
% img_name_set = {'flower1.bmp', 'flower2.bmp'};% image name
% img_name_set = {'62_1_blurred.png','62_2_blurred.png'};

img_name = 'flower1.bmp';

[B_org,B] = preparedata(img_name);
%B = get_levin_data();
ksz = [35,35];% kernel size for each image
%% estimate kernel
%input: B in [0,1]
opts.lambda = 2e-2;% weight for L0
opts.min_lambda = 0;
opts.xk_iter = 30;
opts.iter_step = 1;

opts.verbose = 1;% visualize 
opts.k_threshold = 1/20;
opts.gamma = 1;% weight for kernel
opts.gamma_correct = 1;
%%
if opts.gamma_correct ~= 1
    B = B.^opts.gamma_correct;
end

lambda_scale = initial_weights(B);
opts.lambda = opts.lambda*lambda_scale; 
opts.min_lambda = opts.min_lambda*lambda_scale;


K = kernel_initialize_mul(ksz);

[K,L,opts] = estimate_kernel_single(B,K,opts); 

%center the kernel
K = adjust_psf_center_mul(K);
K = kernel_normalize_mul(K);  
% set elements below threshold to 0
if opts.k_threshold > 0
    K = kernel_normalize_mul(K,opts.k_threshold);
end
figure(1),imshow(K);
figure(2),imshow(L); 
title('Final result');
%% final deblur


