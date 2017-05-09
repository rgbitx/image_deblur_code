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

thresholdR = 1;
thresholdS = 2;


% B=imgaussfilt(B,1.6);

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
  tic
  for itr=maxitr+1:-1:1
      cret=retv(itr);
      Bp=downSmpImC(B,cret);
      if itr == maxitr+1
          I = Bp;
      else
          I = imresize(I, size(Bp) , 'bicubic');%%%%%%%
      end
      
%       Bp=imgaussfilt(Bp,1.2);
      
      r = gradient_confidence_full(Bp,window_size);
      
      Bx = conv2(Bp, [-1,1;0,0], 'valid'); %vertical edges
      By = conv2(Bp, [-1,0;1,0], 'valid'); 
      
      %For each level
      i=1;
      while i<=5
          % evolution
                    
%           r= exp(-r.^(0.8));        
%           structure = structure_adaptive_map( I,lambda_texture*r, 100);          
%           I_= shock(structure,iter,dt,h,'org');
%           structure = I_;
          
          structure= shock(I,iter,dt,h,'org');
                 
          %%%%%%%%%%%compute the gradient of I
          
          I_x = conv2(structure, [-1,1;0,0], 'valid'); %vertical edges
          I_y = conv2(structure, [-1,0;1,0], 'valid'); % horizontal edges         
          I_mag = sqrt(I_x.^2+I_y.^2);
          
%           I_x = I_x.*heaviside_function(I_mag,edge_thresh);
%           I_y = I_y.*heaviside_function(I_mag,edge_thresh);
          
          [r1,c1]=size(I_x);
          MaskR = heaviside_function(r(1:r1,1:c1), thresholdR);           
          I_x = I_x.*heaviside_function(MaskR.*I_mag,thresholdS);
          I_y = I_y.*heaviside_function(MaskR.*I_mag,thresholdS);
          
          %exclude the edge
          I_x(1:2,:) =0;I_x(end-1:end,:) = 0;I_x(:,1:2) =0;I_x(:,end-1:end) = 0;
          I_y(1:2,:) =0;I_y(end-1:end,:) = 0;I_y(:,1:2) =0;I_y(:,end-1:end) = 0;
                       
          if (any(I_x(:))==0) || (any(I_y(:))==0)
              thresholdR = thresholdR/1.1;
              thresholdS = thresholdS/1.1;
              continue;
          end
          
          
          ks=estimate_psf(Bx, By, I_x, I_y, 2, size(kernel));
          
         %% center the kernel
          ks = adjust_psf_center(ks);
          ks(ks(:)<0) = 0;
          ks = ks./sum(ks(:));
          kernel=ks;         
          
          
          lambda_pixel = 4e-3;
          lambda_grad = 4e-4;
%           I = deconv_ansio_L1(Bp,kernel,lambda_smooth,100);
%           I = L0Deblur_whole(Bp, kernel, lambda_pixel, lambda_grad, 2.0);

          I = L0Restoration(Bp, kernel, lambda_grad, 2.0);
          
          I(I<0) = 0;
          I(I>1) = 1;
          
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
          thresholdR = thresholdR/1.1;
          thresholdS = thresholdS/1.1;
          i=i+1;
      end
      
      if (itr>1)
          kernel=resizeKer(kernel,1/ret,k1list(itr-1),k2list(itr-1));
%           kernel = adjust_psf_center(kernel);
%           kernel(kernel<0) = 0;
%           kernel = kernel./sum(kernel(:));
      end
      
  end
  toc
kernel = kernel/sum(kernel(:));
latent = [];

%% kernel refine

%%%%get epsilon_s
% threshold for epsilon_s
[width,height] = size(kernel);
k = kernel;

% threshD = norm(kernel,inf)/(2*width*it);

for it=1:5
    
threshD = norm(k,inf)/(2*width*it);

sorted_k = sort(k(:));
nonzero_index = find(sorted_k>0);
first_nonzero_index = nonzero_index(1);

diff_sorted_k = diff(sorted_k);
all_index = find(diff_sorted_k(first_nonzero_index:end)>threshD);

if isempty(all_index)
   continue; 
end

epsilon_s = sorted_k(all_index(1)+first_nonzero_index);

mask_s = heaviside_function(k,epsilon_s);
mask_s_ = ~mask_s;

s = mask_s .* k;
s_ = mask_s_ .* k;

% set up options for the kernel estimation

k_prev = k;
figure,imagesc(k);
figure,imagesc(s);
figure,imagesc(s_);

k = psf_refine_irls(Bx, By, I_x, I_y, k, s_, 2, size(kernel));
k = adjust_psf_center(k);
k(k(:)<0) = 0;
k = k./sum(k(:));
epsilon_k = norm(k-k_prev)/norm(k_prev)


end

kernel = k;
%% final deblur
fprintf('image deblurring...\n');
for c= 1:size(ori_B,3)
    latent(:,:,c) = deconvSps(ori_B(:,:,c),kernel,0.003,100);
    %% using the edge preserveing model
     %latent(:,:,c) = deconvSps_adaptive_L1(ori_B(:,:,c),kernel,0.003,200,structure);
end


