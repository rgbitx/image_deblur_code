clc;
clear all;
close all;

%2016/06/01 %%%%%%%%%%%%%%%%%%%%%%% Happy Children's Day!!!
% single image  uniform deblurring

% addpath(genpath('cho_code'));
opts.prescale = 1; %%downsampling
opts.xk_iter = 5; %% the iterations
opts.gamma_correct = 1;
opts.k_thresh = 20;
opts.struct = 0.001;
opts.a = 0.01; % 0.01
opts.b = 1.5;  % 1.5 

fn = 'image3.png'; opts.kernel_size = [25 25]; mu = 1500; % 0.1  mu=1000; % dt=0.25

%% use Hyper-Laplacian priors to estimate L

y = imread(fn);
[hei wid dimen] = size(y);
for k = 1:dimen
    y1(:, :, k) = imresize(y(:, :, k), opts.prescale, 'bilinear');
end;
y = y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% % 2014/12/23 修改

% if choose a patch for deblurring by users
isselect = 0; 
if isselect ==1
    figure, imshow(y);
    %tips = msgbox('Please choose the area for deblurring:');
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
[kernel, interim_latent] = deblurring(yg, opts);
toc
y1 = im2double(y1);

%% Final Deblurring:
b = zeros(opts.kernel_size);
bhs = floor(size(b, 1)/2);
for j = 1:size(y,3)
    pad(:,:,j) = padarray(y1(:,:,j), [1 1] * bhs, 'replicate', 'both');
end

% for a = 1:4
%     pad = edgetaper(pad, kernel);
% end

%%%%%%%%%%%%%%%%%%%%%%%%
% pad = shock_filter(pad, 0.3, 1);  % 2014/11/02 修改
%%%%%%%%%%%%%%%%%%%%%%%%

%% Hyper-laplacian deconvolution
for c = 1:size(y,3)
%    
Latent1(:,:,c) = deblurring_adm_aniso_1(pad(:,:,c), kernel,  0.001, 0.8, 0, max(size(kernel(1,2))));  % 去掉了lambda_3，即S的约束项

end
Latent1 = Latent1(bhs + 1 : end - bhs, bhs + 1 : end - bhs, :);

%% TV deconvolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.rho_r   = 2;
opts.beta    = [1 1 0];
opts.print   = true;
opts.alpha   = 0.7;
opts.method  = 'l2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_image = deconvtv(pad, kernel, mu, opts);
Latent2 = L_image.f;
Latent2 = Latent2(bhs + 1 : end - bhs, bhs + 1 : end - bhs, :);

%% the mean value of the two methods
Latent = (Latent1 + Latent2)/2;
figure(2); imshow(Latent);
figure(4),subplot(1,2,1),imshow(Latent1);
subplot(1,2,2),imshow(Latent2);

%% save images
k = kernel - min(kernel(:));
k = k./max(k(:));
imwrite(Latent1,['results\' fn(7:end-4) '_latent1_result.png']);
imwrite(k,['results\' fn(7:end-4) '_kernel.png']);
imwrite(Latent,['results\' fn(7:end-4) '_result_ours.png']);
imwrite(Latent2,['results\' fn(7:end-4) '_latent2_result.png']);

%%

