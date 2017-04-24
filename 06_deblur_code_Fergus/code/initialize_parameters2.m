function [x1,x2] = initialize_parameters2(obs,blur,im,true_blur,true_im,pres,prior,prior_num,mode_im,mode_blur,obs_im,big_blur,spatial_mask,priors,FFT_MODE,COLOR,nLayers)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology


  %%%%%%%%%%%%%
  %%%% Blur
  if strcmp(mode_blur,'direct')
    me2 = blur;
  elseif strcmp(mode_blur,'true')
    me2 = true_blur;
  elseif strcmp(mode_blur,'updown')
    [K,L] = size(true_blur);
    me2 = imresize(imresize(true_blur,0.5,'bilinear'),[K L],'bilinear');
  elseif strcmp(mode_blur,'delta')
    [K,L] = size(true_blur);
    me2 = delta_kernel(K);
  elseif strcmp(mode_blur,'hbar')
    [K,L] = size(true_blur);
    hK = floor(K/2)+1;
    hL = floor(L/2)+1;
    me2 = zeros(K,L);
    me2(hK,hL) = 1;
    me2(hK,hL-1) = 1;
    me2(hK,hL+1) = 1;

  elseif strcmp(mode_blur,'vbar')
    [K,L] = size(true_blur);
    hK = floor(K/2)+1;
    hL = floor(L/2)+1;
    me2 = zeros(K,L);
    me2(hK,hL) = 1;
    me2(hK-1,hL) = 1;
    me2(hK+1,hL) = 1;

  elseif strcmp(mode_blur,'star')
    [K,L] = size(true_blur);
    hK = floor(K/2)+1;
    hL = floor(L/2)+1;
    me2 = zeros(K,L);
    me2(hK-1,hL+1) = 1;
    me2(hK-1,hL-1) = 1;
    me2(hK+1,hL-1) = 1;
    me2(hK+1,hL+1) = 1;


  elseif strcmp(mode_blur,'random')
    [M,N] = size(true_blur);
    me2 = rand(M,N);
  elseif strcmp(mode_blur,'variational')
    me2 = blur;
    %%TODO
  else
    error foo
  end
  
  
  %%%%%%%%%%%%%%%
  %%%% Image

  if COLOR
    C = 3;
  else
    C = 1;
  end

  spatial_mask = spatial_mask(:);
  
  if strcmp(mode_im,'direct')
    mx2 = im;
  elseif strcmp(mode_im,'true')
    mx2 = true_im;
  elseif strcmp(mode_im,'slight_blur_obs')

    [M,N] = size(obs);
    f = [1 2 1; 2 4 2; 1 2 1]/16;

    for c=1:C
      mx2(:,:,c) = real(ifft2(fft2(f,M,N).*fft2(obs(:,:,c),M,N)));
    end
    
  elseif strcmp(mode_im,'lucy')

    obs2 = edgetaper(obs_im,big_blur);
    im_tmp = deconvlucy(obs2,big_blur);
    
    for c=1:C
      im_x(:,:,c) = conv2(im_tmp(:,:,c),[1 -1],'valid');
      im_y(:,:,c) = conv2(im_tmp(:,:,c),[1 -1]','valid');
    end
    
    [M,N] = size(obs);
    
    im_xs = imresize(im_x,[M N/2]+2,'bilinear');
    im_ys = imresize(im_y,[M N/2]+2,'bilinear');
   
    mx2 = zeros(M,N,C);
    mx2(1:M,1:N/2,:) = im_xs(2:end-1,2:end-1,:);
    mx2(1:M,N/2+1:N,:) = im_ys(2:end-1,2:end-1,:);
    
  elseif strcmp(mode_im,'reg')
    obs2 = edgetaper(obs,blur);
    mx2 = deconvreg(obs2,blur);
  elseif strcmp(mode_im,'lucy_true')
    obs2 = edgetaper(obs,true_blur);
    mx2 = deconvlucy(obs2,true_blur,20);
  elseif strcmp(mode_im,'reg_true')
    obs2 = edgetaper(obs,true_blur);
    mx2 = deconvreg(obs2,true_blur);
  elseif strcmp(mode_im,'updown')
    [M,N] = size(true_im);
    mx2 = imresize(imresize(true_im,0.5,'bilinear'),[M N],'bilinear');
  elseif strcmp(mode_im,'random')
    %%% generate image using laplacian distribution of 
    SCALE_PARAMETER = [7 6 4];
     [M,N] = size(im);
  
    tmp1 = rand(M,N,C);
    tmp2 = rand(M,N,C);
    
    mx2 = SCALE_PARAMETER(3) * (-log(tmp1) + log(tmp2));
    
  elseif strcmp(mode_im,'variational')
    MAX_ITERATIONS = 5000;
    [M,N] = size(im);
    [K,L] = size(blur);
    %norm_blur = blur / sum(blur(:));
    norm_blur = me2 / sum(me2(:));

    %%% n.b. for SIGGRAPH only 2 blur components were used.
    %%% set dimensions(2,6)=1 to get back to SSG & IMA runs
    
    dimensions = [1 1   1  0  0  1;
                  1 K*L 4  1  0  0;
                  C M*N 4  0  1  1];
 
    
    if FFT_MODE
      
      I = M*2; J = N*2;  
      Dpf = zeros(I,J,C);
      Dpf(K:M,L:N/2,:) = 1;
      Dpf(K:M,L+N/2:N,:) = 1; %% y-plane
      
      Df = padarray(obs,[M N],0,'post');
      
   else
      
      I = M; J = N;
      hK = floor(K/2); hL = floor(L/2);
      
      Dpf = zeros(I,J,C);
      Dpf(hK+1:M-hK,hL+1:N/2-hL,:) = 1;
      Dpf(hK+1:M-hK,N/2+hL+1:N-hL,:) = 1;
      
      shift_kernel = zeros(K,L);
      shift_kernel(1,1)=1;
    
      for c=1:C
        Df(:,:,c) = conv2(obs(:,:,c),shift_kernel,'same');
      end
    end
    
    pres_vector = ones(1,length(norm_blur(:))+length(im(:))+1) * pres;
    q=find(spatial_mask(:));
    pres_vector(q+1+length(norm_blur(:))) = spatial_mask(q)';
            
    dummy_blur_mask = zeros(dimensions(2,3),length(norm_blur(:)));
    
    %%% Make vectors
    xx1 = [0 norm_blur(:)' im(:)'] .* pres_vector;  
    xx2 = pres_vector;

 
    [ensemble,D_log,gamma_log]=train_ensemble_main6(dimensions,xx1,xx2,'','',[1e-4 0 1 0 0 MAX_ITERATIONS 0],Df,Dpf,I,J,K,L,M,N,priors,FFT_MODE,dummy_blur_mask,[1-(spatial_mask>0)]);    
     
    mx2 = reshape(train_ensemble_get(3,dimensions,ensemble.mx),M,N,C);

  elseif strcmp(mode_im,'greenspan')  
  
      %%% Use default parameters
      s = create_greenspan_settings;
      s.factor = 0;

      [en,lo] = greenspan(obs,s);
      mx2 = en;
      %%% get size of estimated image at current scale
      %[M,N] = size(im);

      %%% downsample obs_im (must work in intensity space)
      %obs_small = imresize(obs_im,([M N/2]/2)+1,'bilinear');
      
      %%% use greenspan
      %[en,lo] = greenspan(obs_small,s);
      
      %%% Now go back to gradient space
      %x_plane = conv2(en,[1 -1],'valid');
      %y_plane = conv2(en,[1 -1]','valid');
      
      %mx2 = [x_plane(1:M,1:N/2), y_plane(1:M,1:N/2)];

  else
    error foo
  end
  


  %%% ensure blur is normalized to 1
  me2 = me2 / sum(me2(:));
%  mx2  = mx2 * std(obs(:))/std(mx2(:));


  %%% do spatial masking
  pres_vector = ones(1,length(me2(:))+length(mx2(:))+nLayers) * pres;
  q=find(spatial_mask(:));
  pres_vector(q+1+length(me2(:))) = spatial_mask(q)';
 
  %%% Make vectors
  x1 = [zeros(1,nLayers) me2(:)' mx2(:)'] .* pres_vector;  
  x2 = pres_vector;
