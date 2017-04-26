function [k, S] = sub_deblurring(blur_B, k, threshold, opts)

%2016/5/31 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Input:
% @blur_B: input blurred image
% @k: blur kernel

% Ouput:
% @k: estimated blur kernel
% @S: intermediate latent image
%

% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = size(blur_B,1);    W = size(blur_B,2);
[k1,k2]= size(k);
k_max = max(k1, k2);
bhs = floor(k_max/2);
blur_B_tmp = blur_B(1:H,1:W,:);
Bx = conv2(blur_B_tmp, dx, 'valid');
By = conv2(blur_B_tmp, dy, 'valid');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = blur_B;

for iter = 1:opts.xk_iter
    %% latent image deblurring
    % salient structure extraction    
    blur_B = padarray(blur_B, [1 1] * bhs, 'replicate', 'both');    
    S = deblurring_adm_aniso_1(blur_B, k,  0.001, 0.5, 0.05, max(size(k(1,2))));
    S = S(bhs + 1 : end - bhs, bhs + 1 : end - bhs, :);
    blur_B = blur_B(bhs + 1 : end - bhs, bhs + 1 : end - bhs, :);
    blur_B = blur_B(1:H,1:W,:);
    S = wrap_boundary_liu(S, opt_fft_size([H W]+size(k)-1));
    S = S(1:H,1:W,:);
%     figure(5), imshow(S);
    S = tsmooth(S, opts.a, opts.b);
%     figure(6),imshow(S);
    opts.a = opts.a*1.1;
    opts.b = opts.b*1.1;
    dt=1; h=1;
    iter_shock=10;  % 100
    S = shock(S,iter_shock,dt,h,'org');
    [latent_x, latent_y, threshold]= threshold_pxpy_v1(S,max(size(k)),threshold);
    k_prev = k;
    
    % kernel estimation
    tic
    k = estimate_psf(Bx, By, latent_x, latent_y, 2, size(k_prev));
    toc
   %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Hu et al
    CC = bwconncomp(k,8);
    for ii=1:CC.NumObjects
        currsum=sum(k(CC.PixelIdxList{ii}));
        if currsum<.1
            k(CC.PixelIdxList{ii}) = 0;
        end
    end
    k(k<0) = 0;
    k=k/sum(k(:));

    figure(1);
    S(S<0) = 0;
    S(S>1) = 1;
    subplot(1,3,1); imshow(blur_B,[]); title('Blurred image');
    subplot(1,3,2); imshow(S,[]);title('Interim latent image');
    subplot(1,3,3); imshow(k,[]);title('Estimated kernel');
    end;
k(k<0) = 0;
k = k ./ sum(k(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
