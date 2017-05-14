clear,clc,close all;
addpath('cho_code');
%% get input
% img_name_set = {'book1.png', 'book2.png','book3.png'};
% img_name_set = {'stone1.png', 'stone2.png'};
img_name_set = {'flower1.bmp', 'flower2.bmp'};% image name
% img_name_set = {'62_1_blurred.png','62_2_blurred.png'};
[B_org,B] = preparedata(img_name_set);
% B = get_levin_data();
ksz = [21,21];% kernel size for each image
%% estimate kernel
% for levin data: lambda = 1e-3
opts.lambda = 5e-2;% weight for L0
opts.min_lambda = 1e-4;
opts.xk_iter = 20;
opts.iter_step = 5;

opts.verbose = 0;% visualize 
opts.min_ksz = 3;
opts.scale_factor = sqrt(2);
opts.k_threshold = 1/20;
opts.gamma = 1;% weight for kernel
opts.gamma_correct = 1;
%%
if opts.gamma_correct ~= 1
    B = B.^opts.gamma_correct;
end
Ks = kernel_initialize_mul(ksz,size(B,3));
[nLevels,mbsz,mksz] = prepare_multi_scale([size(B,1),size(B,2)],ksz,opts.min_ksz,...
                                           opts.scale_factor);
for iLevel = 1:nLevels    
    Ks = imresize(Ks,mksz(iLevel,:),'bilinear');
    Ks = kernel_normalize_mul(Ks);    
    Bs = imresize(B,mbsz(iLevel,:),'bilinear');    
    %% Initialize the parameter
    % update kernel
    [Ks,Ls,opts] = estimate_kernel_single(Bs,Ks,opts);
    figure(1),subplot(1,2,1),imshow(mat2gray(Ks(:,:,1))),
    subplot(1,2,2),imshow(mat2gray(Ks(:,:,2)));
    figure(2),imshow(mat2gray(Ls));
    title(sprintf('level=%d',iLevel));
    pause(0.01);    
    %center the kernel
    %Ks = adjust_psf_center_mul(Ks);
    %Ks = kernel_normalize_mul(Ks);  
end%iLevel
% set elements below threshold to 0
K = Ks;
if opts.k_threshold > 0
    K = kernel_normalize_mul(K,opts.k_threshold);
end
%% final deblur