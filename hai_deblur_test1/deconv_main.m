function [kernel, latent] = deconv_main(B, ori_B, lambda_kernel, lambda_smooth,...
    lambda_texture, window_size, k1, k2, lambda_kernel_smooth, display)
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
 %%initialize the kernel 
 kernel = zeros(k1,k2);
maxitr=max(floor(log(5/min(k1,k2))/log(ret)),0);

retv=ret.^[0:maxitr];

k1list=ceil(k1*retv);
k1list=k1list+(mod(k1list,2)==0);
k2list=ceil(k2*retv);
k2list=k2list+(mod(k2list,2)==0);

cret=retv(end);
kernel=resizeKer(kernel,cret,k1list(end),k2list(end));
% parameters
dt=1; h=1;
iter=100;      % iterations
% %fixed the size of kernel, result seems better, but time consuming
% %   kernel = zeros(k1,k2);
  I_X = conv2(B, [-1,1;0,0], 'valid'); %vertical edges
  I_Y = conv2(B, [-1,0;1,0], 'valid'); %
  I_mag = sqrt(I_X.^2+I_Y.^2);
  [k1,k2]= size(kernel);
  edge_thresh = determine_truck(I_X, I_Y, I_mag,k1,k2);
  clear I_X I_Y I_mag;
  for itr=maxitr+1:-1:1
      cret=retv(itr);
      Bp=downSmpImC(B,cret);
      if itr == maxitr+1
          I = Bp;
      else
          I = imresize(I, size(Bp) , 'bicubic');%%%%%%%
      end
  %For each level
      for i = 1:5
          % evolution
          
          r = gradient_confidence_full(Bp,window_size);
          r= exp(-r.^(0.8));
          
          structure = structure_adaptive_map( I,lambda_texture*r, 100);
          I_= shock(structure,iter,dt,h,'org');
          structure = I_;
          
          %%%%%%%%%%%compute the gradient of I
          I_x = conv2(structure, [-1,1;0,0], 'valid'); %vertical edges
          I_y = conv2(structure, [-1,0;1,0], 'valid'); % horizontal edges
          Bx = conv2(Bp, [-1,1;0,0], 'valid'); %vertical edges
          By = conv2(Bp, [-1,0;1,0], 'valid'); 
          I_mag = sqrt(I_x.^2+I_y.^2);
          I_x = I_x.*heaviside_function(I_mag,edge_thresh);
          I_y = I_y.*heaviside_function(I_mag,edge_thresh);
          %exclude the edge
          I_x(1:2,:) =0;I_x(end-1:end,:) = 0;I_x(:,1:2) =0;I_x(:,end-1:end) = 0;
          I_y(1:2,:) =0;I_y(end-1:end,:) = 0;I_y(:,1:2) =0;I_y(:,end-1:end) = 0;
                       
          tic
          kernel=estimate_psf(Bx, By, I_x, I_y, 2, size(kernel));
          toc
          
          tic
          I = deconv_ansio_L1(Bp,kernel,lambda_smooth,100);
          toc
          
          fprintf('%d iterations of pyramid %d is done\n', i, itr);
          %% for test
          if display==1
              figure(1);subplot(221); imshow(Bp,[]);title('Input')
              subplot(222);imshow(structure,[]);title('structure')
              subplot(223); imshow(I,[]); title('deblurring image');
              subplot(224); imshow(fliplr(flipud(kernel)),[]); title('kernel');
          end
          dt = dt/1.1;
          lambda_texture = lambda_texture/1.1;
          edge_thresh = edge_thresh/1.1;
      end
      if (itr>1)
          kernel=resizeKer(kernel,1/ret,k1list(itr-1),k2list(itr-1));
          kernel = adjust_psf_center(kernel);
          kernel(kernel<0) = 0;
          kernel = kernel./sum(kernel(:));
      end
  end
kernel = kernel/sum(kernel(:));
latent = [];
fprintf('image deblurring...\n');
for c= 1:size(ori_B,3)
    latent(:,:,c) = deconvSps(ori_B(:,:,c),kernel,0.003,100);
    %% using the edge preserveing model
     %latent(:,:,c) = deconvSps_adaptive_L1(ori_B(:,:,c),kernel,0.003,200,structure);
end


