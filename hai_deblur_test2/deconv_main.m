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

lambda_grad = 4e-4;

%%initialize the kernel 
k1 = opt.kernel_sizeh;
k2 = opt.kernel_sizew;
kernel = zeros(k1,k2);
maxitr=max(floor(log(5/min(k1,k2))/log(ret)),0);

thresholdR = 0.5;
thresholdS = 2*thresholdR;

% B=imgaussfilt(B,1.6);

retv=ret.^[0:maxitr];

k1list=ceil(k1*retv);
k1list=k1list+(mod(k1list,2)==0);
k2list=ceil(k2*retv);
k2list=k2list+(mod(k2list,2)==0);

cret=retv(end);
kernel=resizeKer(kernel,cret,k1list(end),k2list(end));

% I_X = conv2(B, [-1,1;0,0], 'valid'); %vertical edges
% I_Y = conv2(B, [-1,0;1,0], 'valid'); %
% I_mag = sqrt(I_X.^2+I_Y.^2);
% [k1,k2]= size(kernel);
% edge_thresh = determine_truck(I_X, I_Y, I_mag,k1,k2);
% clear I_X I_Y I_mag;
% parameters
dt=1; h=1;
iter=100;      % iterations
% get proper thresholdD and thresholdS
maxCnt = 200;
  
  for itr=maxitr+1:-1:1
              
      cret=retv(itr);
      Bp=downSmpImC(B,cret);
      if itr == maxitr+1
          I = Bp;
      else
          I = imresize(I, size(Bp) , 'bicubic');%%%%%%%
      end
   
      r = gradient_confidence_full(Bp,opt.window_size);
      
      Bx = conv2(Bp, [-1,1;0,0], 'valid'); %vertical edges
      By = conv2(Bp, [-1,0;1,0], 'valid'); 
            
      %For each level
      i=1;
      while i<=5
          % evolution
                                    
          structure= shock(I,iter,dt,h,'org');
                 
          %%%%%%%%%%%compute the gradient of I
          
          I_x0 = conv2(structure, [-1,1;0,0], 'valid'); %vertical edges
          I_y0 = conv2(structure, [-1,0;1,0], 'valid'); % horizontal edges         
          I_mag = sqrt(I_x0.^2+I_y0.^2);
          
          
          [r1,c1]=size(I_x0);
          MaskR = heaviside_function(r(1:r1,1:c1), thresholdR);           
          I_x = I_x0.*heaviside_function(MaskR.*I_mag,thresholdS);
          I_y = I_y0.*heaviside_function(MaskR.*I_mag,thresholdS);
          
          %exclude the edge
          I_x(1:2,:) =0;I_x(end-1:end,:) = 0;I_x(:,1:2) =0;I_x(:,end-1:end) = 0;
          I_y(1:2,:) =0;I_y(end-1:end,:) = 0;I_y(:,1:2) =0;I_y(:,end-1:end) = 0;
          
          % make I_x and I_y are nonzero
          cnt = 1;
          while cnt < maxCnt && ((any(I_x(:))==0) || (any(I_y(:))==0))
              
              thresholdR = thresholdR/1.1;
              thresholdS = thresholdS/1.1;

              MaskR = heaviside_function(r(1:r1,1:c1), thresholdR);           
              I_x = I_x0.*heaviside_function(MaskR.*I_mag,thresholdS);
              I_y = I_y0.*heaviside_function(MaskR.*I_mag,thresholdS);

              cnt = cnt + 1;
          end
          if cnt==maxCnt
              break;
          end
          
                  
          kernel=estimate_psf(Bx, By, I_x, I_y, 20, size(kernel));
          
         %% center the kernel
          kernel = adjust_psf_center(kernel);
          kernel(kernel(:)<0) = 0;
          kernel = kernel./sum(kernel(:));
          
          
          I = L0Restoration(Bp, kernel, lambda_grad);
          
          I(I<0) = 0;
          I(I>1) = 1;
                                 
          fprintf('%d iterations of pyramid %d is done\n', i, itr);
          %% for test
          if opt.display==1
              figure(1);subplot(221); imshow(Bp,[]);title('Input')
              subplot(222);imshow(structure,[]);title('structure')
              subplot(223); imshow(I,[]); title('deblurring image');
              subplot(224); imshow(fliplr(flipud(kernel)),[]); title('kernel');
          end
          dt = dt/1.1;
%           lambda_grad = max(lambda_grad/1.1,4e-4);
         
          thresholdR = thresholdR/1.1;
          thresholdS = thresholdS/1.1;
          i=i+1;
      end
      
      if (itr>1)
          kernel=resizeKer(kernel,1/ret,k1list(itr-1),k2list(itr-1));
          kernel = adjust_psf_center(kernel);
          kernel(kernel<0) = 0;
          kernel = kernel./sum(kernel(:));
      end
      
  end
thresholdR
thresholdS
kernel = kernel/sum(kernel(:));
latent = [];

%% final deblur
fprintf('image deblurring...\n');
opts.nb_lambda = 3000;
opts.nb_alpha = 1.0;
bhs = floor(size(kernel,1)/2);

for c= 1:size(ori_B,3)
    ypad = padarray(ori_B(:, :, c), [1 1] * bhs, 'replicate', 'both');
    tmp = fast_deconv_bregman(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
    latent(:, :, c) = tmp(bhs + 1 : end - bhs, bhs + 1 : end - bhs);
%     latent(:,:,c) = deconvSps(ori_B(:,:,c),kernel,0.003,100);
    %% using the edge preserveing model
%      latent(:,:,c) = deconvSps_adaptive_L1(ori_B(:,:,c),kernel,0.003,200,structure);
end



