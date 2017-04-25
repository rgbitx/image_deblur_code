function [I, ex, ssde, kernel] = deconv_main_multiple(B, ori_B, lambda_kernel, lambda_smooth,...
    lambda_texture, window_size, k1, k2, lambda_kernel_smooth, display, x)
%Created by Jinshan Pan
%More details can be found:
%Kernel Estimation from Salient Structure for Robust Motion Deblurring
ret = sqrt(0.5);
 %%initialize the kernel 
 kernel = zeros(k1,k2);
maxitr=max(floor(log(5/min(k1,k2))/log(ret)),0);
%maxitr = 1;
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
      I_Y = conv2(B, [-1,0;1,0], 'valid'); %vertical edges
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
  %every level solver
      for i = 1:5
          % evolution
          r = gradient_confidence_full(Bp,window_size);
          r= exp(-r.^(0.8));
          structure = structure_adaptive_map( I,lambda_texture*r, 100,0.95);
          %structure  = I;
          I_= shock(structure,iter,dt,h,'org');
          structure = I_;
% %           %%for test
% %           figure(1);subplot(221); imshow(Bp,[]);title('original')
% %           subplot(222);imshow(structure,[]);title('structure')

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
          if size(I_x) ==size(Bx)
              [N1,N2]=size(I_x);
              hfs1_x1=floor((size(kernel,2)-1)/2);
              hfs1_x2=ceil((size(kernel,2)-1)/2);
              hfs1_y1=floor((size(kernel,1)-1)/2);
              hfs1_y2=ceil((size(kernel,1)-1)/2);
              hfs_x1=hfs1_x1;
              hfs_x2=hfs1_x2;
              hfs_y1=hfs1_y1;
              hfs_y2=hfs1_y2;
              N2=N2+hfs_x1+hfs_x2;
              N1=N1+hfs_y1+hfs_y2;
              N=N1*N2;
              mask=zeros(N1,N2);
              mask(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=1;
              tI_x=I_x;
              I_x=zeros(N1,N2);
              I_x(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=tI_x;
              tI_y=I_y;
              I_y=zeros(N1,N2);
              I_y(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=tI_y;
          end
          [k1,k2] = size(kernel);
          A=zeros(k1*k2);
          b=zeros(k1*k2,1);
          for j=1:size(I,3)
              A = A+ getAutoCor(I_x(:,:,j),k1,k2);
              b = b+getCory(I_x(:,:,j),Bx(:,:,j),k1,k2);
              A = A+ getAutoCor(I_y(:,:,j),k1,k2);
              b = b+getCory(I_y(:,:,j),By(:,:,j),k1,k2);
          end
          %% will be replaced by our fast algorithm 
          kernel =solve_kernel_irls(A,b,k1,k2,lambda_kernel, lambda_kernel_smooth);
          %%
          %ks = size(kernel)-1;
          kernel  = Cho_correct(kernel);
          %kernel = kernel(ks+1:end-ks, ks+1:end-ks);
          kernel(kernel<0) = 0;
          kernel = kernel/sum(kernel(:));
          %% will be replaced by our fast algorithm 
          %I = deconv_ansio_L1(Bp,(kernel),lambda_smooth,30);
          %I = deblurring_lr_whole(Bp, kernel, 0.01);
          %I = deblurring_lr_whole_with_edge(Bp, kernel, 0.01, structure, edge_thresh);
          I = deblurring_lr_whole_with_edge_without_lr(Bp, kernel, 0.01, structure, edge_thresh);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %Lin Zhouchen's suggestion
          fprintf('%d iterations of pyramid %d is done\n', i, itr);
          %% for test
%           if display==1
%               figure(1);subplot(221); imshow(Bp,[]);title('original')
%               subplot(222);imshow(structure,[]);title('structure')
%               subplot(223); imshow(I,[]); title('deblurring image');
%               subplot(224); imshow(fliplr(flipud(kernel)),[]); title('kernel');
%           end
          dt = dt/1.1;
          lambda_texture = lambda_texture/1.1;
          edge_thresh = edge_thresh/1.1;
      end
      if (itr>1)
          kernel=resizeKer(kernel,1/ret,k1list(itr-1),k2list(itr-1));
      end
  end
  
%ori_B = double(ori_B)/255;
kernel = kernel/sum(kernel(:));
fprintf('image deblurring...\n');
[ex]=deconvSps(ori_B,kernel,0.0068,70);
% ex = 1;
[ssde]=comp_upto_shift(ex,x);
% ssde = 1;

