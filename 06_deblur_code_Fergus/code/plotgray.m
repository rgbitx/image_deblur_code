function plotgray(tmp_name)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology
  
  close all;
  EXPORT_DIR = './';
  % plot initialisation of inference as well as the output of each
  % inference stage
  PLOT_INIT = 0;
  
  load([tmp_name,'.mat']);

  if exist(['tmp_',tmp_name,'.mat'])
    load(['tmp_',tmp_name]);
  end
 
  if (SYNTHETIC==1)
    fprintf('WARNING - Synthetic image\n');
  else
    fprintf('Real image\n');
  end
 
  if PLOT_INIT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Initial image
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure; 
  for a=1:length(mx_init)
    ultimateSubplot(NUM_SCALES,1,a,[],0.1);
    dx = mx_init{a}(:,1:size(mx_init{a},2)/2);
    dy = mx_init{a}(:,size(mx_init{a},2)/2+1:end);
    dx(:,1) = 0; dx(1,:) = 0; dx(:,end) = 0; dx(end,:) = 0;
    dy(:,1) = 0; dy(1,:) = 0; dy(:,end) = 0; dy(end,:) = 0;    
    [obs_tmp,k] = reconsEdge3(dx,dy); 
    imagesc(obs_tmp); colormap(gray); colorbar;
    title(['Init - Scale ',num2str(a)]);
  end 
  set(h,'Position',[ 66         746        1452         354]);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Estimated image
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure; 
  for a=1:length(mx_est)
    ultimateSubplot(NUM_SCALES,1,a,[],0.1);
    dx = mx_est{a}(:,1:size(mx_est{a},2)/2);
    dy = mx_est{a}(:,size(mx_est{a},2)/2+1:end);
    dx(:,1) = 0; dx(1,:) = 0; dx(:,end) = 0; dx(end,:) = 0;
    dy(:,1) = 0; dy(1,:) = 0; dy(:,end) = 0; dy(end,:) = 0;    
    [obs_tmp,k] = reconsEdge3(dx,dy); 
    imagesc(obs_tmp); colormap(gray); colorbar;
    title(['Est - Scale ',num2str(a)]);
  end
  set(h,'Position',[65         402        1451         309]);
  
  if PLOT_INIT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Initial kernels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure; 
  for a=1:length(me_init)
    ultimateSubplot(NUM_SCALES,1,a,[],0.1);
    imagesc(me_init{a}); colormap(gray); colorbar;
    title(['Init - Scale ',num2str(a)]);
  end
  set(h,'Position',[  64         484        1494         237]);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Estimated kernels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure; 
  for a=1:length(me_est)
    ultimateSubplot(NUM_SCALES,1,a,[],0.1);
    imagesc(me_est{a}); colormap(gray); 
    title(['Est - Scale ',num2str(a)]); q=caxis; q(1)=0; caxis(q);colorbar;
  end
  set(h,'Position',[ 63         237        1495         248]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Final patch gallery
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (length(mx_est)==NUM_SCALES)
    h=figure; 
    ultimateSubplot(2,1,1,[],0.1);
    pp = obs_im_all(PATCH_LOCATION(2):PATCH_LOCATION(2)+PATCH_SIZE(2),PATCH_LOCATION(1):PATCH_LOCATION(1)+PATCH_SIZE(1));
    imagesc(pp); colormap(gray); colorbar;
    title('Original blurry patch');
    ultimateSubplot(2,1,2,[],0.1);
    g_obs_tmp = fix_image(obs_tmp,pp);
    imagesc(g_obs_tmp); colormap(gray); colorbar;
    title('Reconstructed patch');
    set(h,'Position',[16         639        1477         420]);
  
    figure; imagesc(obs_im_orig); title('Original blurry image');

    if exist('deblurred_im')
      figure; imagesc(deblurred_im); title('Deblurred image');
    end
    
    if exist('kernel_out') 
      figure; imagesc(kernel_out); axis square; colormap(gray);
      title('Final kernel');
    end
    
  end  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% True blur kernels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if SYNTHETIC
    h=figure; 
    ultimateSubplot(1,1,1,[],0.1);
    imagesc(blur_kernel); colormap(gray); colorbar;
    title(['True kernel']);
  end
  
