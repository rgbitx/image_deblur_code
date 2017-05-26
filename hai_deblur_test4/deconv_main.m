function [kernel, latent] = deconv_main(B, ori_B, opt)
%% 
%%   Parameters: 
%   B: Blurred image (grayscale).
%   ori_B: Original input image.
%   lambda_kernel: Weight for the kernel                     
%   lambda_kernel: Weight for the image     
%   lambda_texture: Weight in equation (3)
%   window_size: The patch size defined in equation (4)
%   k1, k2: The kernel size
%   lambda_kernel_smooth: Weight for l0 norm in equation (10)
%   display: (1/0): display the intermediate results
%
%   The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Risheng Liu, Zhixun Su, and Xianfeng Gu,
%        Kernel estimation from salient structure for robust motion deblurring,
%        Signal Processing: Image Communication, 2013
%  
%   The code and the algorithm are for non-comercial use only.
%   Author: Jinshan Pan (jspan@mail.dlut.edu.cn)
%   Date  : 03/27/2013
%   Copyright 2013, Dalian University of Technology.

ret = sqrt(0.5);


lambda_l0 = 3.8e-4;
lamda_kernel = 5;
%%initialize the kernel 
k1 = opt.kernel_sizeh;
k2 = opt.kernel_sizew;
kernel = zeros(k1,k2);
maxitr=max(floor(log(5/min(k1,k2))/log(ret)),0);

retv=ret.^[0:maxitr];

k1list=ceil(k1*retv);
k1list=k1list+(mod(k1list,2)==0);
k2list=ceil(k2*retv);
k2list=k2list+(mod(k2list,2)==0);

last = maxitr+1;
first = 1;
for itr=last:-1:first
    
    %%%
    if (itr == last)
      %%
      % at coarsest level, initialize kernel
      kernel = init_kernel(k1list(itr));
      kernel = adjust_psf_center(kernel);
%       k1 = k1list(s);
%       k2 = k1; % always square kernel assumed
    else
    % upsample kernel from previous level to next finer level
%     k1 = k1list(s);
%     k2 = k1; % always square kernel assumed
    
    % resize kernel from previous level
      kernel = resizeKer(kernel,1/ret,k1list(itr),k2list(itr));
      kernel = adjust_psf_center(kernel);
      kernel(kernel<0) = 0;
      kernel = kernel./sum(kernel(:));
    end
    
    %%%%  
      cret=retv(itr);
      Bp=downSmpImC(B,cret);
      
      Bx = conv2(Bp, [-1,1;0,0], 'valid'); %vertical edges
      By = conv2(Bp, [-1,0;1,0], 'valid'); 
      
      if itr>=3
        i_max = 3;
      else
        i_max = itr;
      end
      
%       i_max = 1;
      
      %For each level
      i=1;
      while i<=i_max
          % evolution
          I = L0Restoration(Bp,kernel,lambda_l0);
    
          I(I<0) = 0;
          I(I>1) = 1;
                  
          structure= I;
          
          %%%%%%%%%%%compute the gradient of I          
          I_x = conv2(structure, [-1,1;0,0], 'valid'); %vertical edges
          I_y = conv2(structure, [-1,0;1,0], 'valid'); % horizontal edges         
         
          %exclude the edge
          I_x(1:2,:) =0;I_x(end-1:end,:) = 0;I_x(:,1:2) =0;I_x(:,end-1:end) = 0;
          I_y(1:2,:) =0;I_y(end-1:end,:) = 0;I_y(:,1:2) =0;I_y(:,end-1:end) = 0;
                      
          kernel=estimate_psf(Bx, By, I_x, I_y, lamda_kernel, size(kernel));             
                                                
          fprintf('%d iterations of pyramid %d is done\n', i, itr);
          %% for test
%           if opt.display==1
%               figure(1);subplot(221); imshow(Bp,[]);title('Input')
%               subplot(222);imshow(structure,[]);title('structure')
%               subplot(223); imshow(I,[]); title('deblurring image');
%               subplot(224); imshow(rot90(kernel,2),[]); title('kernel');
%               pause(0.001);
%           end
      
          lamda_kernel = lamda_kernel * 1.1;
          
          i=i+1;
      end
      
end



kernel = kernel/sum(kernel(:));
latent = [];

% kernel = imgaussfilt(kernel,1.0);
%% final deblur
fprintf('image deblurring...\n');
% opts.nb_lambda = 3000;
% opts.nb_alpha = 1.0;
opts.nb_lambda = 2000;
opts.nb_alpha = 0.5;
bhs = floor(size(kernel,1)/2);

kernel(kernel==0) = 1e-9;

for c= 1:size(ori_B,3)
    
    ypad = padarray(ori_B(:, :, c), [1 1] * bhs, 'replicate', 'both');   
    for a = 1:1
        ypad = edgetaper(ypad, kernel);
    end
    
%     tmp = fast_deconv_bregman(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
    tmp = fast_deconv(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
    latent(:, :, c) = tmp(bhs + 1 : end - bhs, bhs + 1 : end - bhs);
%     latent(:,:,c) = deconvSps(ori_B(:,:,c),kernel,0.003,100);
    %% using the edge preserveing model
%      latent(:,:,c) = deconvSps_adaptive_L1(ori_B(:,:,c),kernel,0.003,200,structure);
end



end

function [k] = init_kernel(minsize)
  k = zeros(minsize, minsize);
  k((minsize - 1)/2, (minsize - 1)/2:(minsize - 1)/2+1) = 1/2;
end


