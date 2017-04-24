function deblur(CONFIG_FNAME)

% Author: Rob Fergus
% Version: 1.0, distribution code.

% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology
    
%% parse CONFIG_FNAME. If it has an image extension, remove. This assumes
%a .mat file with the same name is present
q = findstr(CONFIG_FNAME,'.');

if isempty(q)
  out_fname = CONFIG_FNAME;
else
  out_fname = CONFIG_FNAME;
  CONFIG_FNAME = CONFIG_FNAME(1:q-1);
end

%%% Loadup all settings
load(char(CONFIG_FNAME));
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULT_GAMMA = 2.2;

if (~exist('SPECIFIC_COLOR_CHANNEL'))
  SPECIFIC_COLOR_CHANNEL = 0;
end

if (~exist('MANUAL_KERNEL'))
  MANUAL_KERNEL = 0;
end

if (~exist('FFT_MODE'))
  FFT_MODE = 1;
end

if (~exist('SATURATION_MASK'))
  SATURATION_MASK = 0;
end

if (~exist('SATURATION_PRESCISION'))
  SATURATION_PRESCISION = 1e1;
end

if (~exist('BLUR_MASK'))
  BLUR_MASK = 0;
end

if (~exist('AUTOMATIC_PATCH'))
  AUTOMATIC_PATCH = 0;
end

if (~exist('BLUR_LOCK'))
  BLUR_LOCK = 0;
end

if (~exist('CONSISTENCY_GRAD_ITS'))
  CONSISTENCY_GRAD_ITS = 0;
end

if (~exist('CAMERA_TYPE'))
  CAMERA_TYPE = 'unknown';
end

if (~exist('LUCY_ITS'))
  LUCY_ITS = 10;
end

if (~exist('KERNEL_THRESHOLD'))
  KERNEL_THRESHOLD = 7;
end

if (~exist('SCALE_OFFSET'))
  SCALE_OFFSET = 0;
end

%%%%% Color not working yet, so turn it off
COLOR = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess observed image 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get size
[obs_imy,obs_imx,obs_imz] = size(obs_im);

%%% Save original image
obs_im_orig = obs_im;

%%% ROB NOTE: shouldn't rgb2gray be AFTER gamma correction?
% Convert to grayscale
if (~COLOR & obs_imz>1)
  %  obs_im = obs_im(:,:,3);
  
  if SPECIFIC_COLOR_CHANNEL
    obs_im = obs_im(:,:,SPECIFIC_COLOR_CHANNEL);
  else
    obs_im = rgb2gray_rob(obs_im);
  end
  
  obs_imz = 1;
end

if (PRESCALE)
  obs_im = imresize(obs_im,PRESCALE,'bilinear');
  obs_im_orig = imresize(obs_im_orig,PRESCALE,'bilinear');
end

% Get size
[obs_imy,obs_imx,obs_imz] = size(obs_im);

% Take gradients
   
  if (GAMMA_CORRECTION~=1)
    %%% gamma correct actual image
    obs_im = (double(obs_im).^(GAMMA_CORRECTION))/(256^(GAMMA_CORRECTION-1));
  else
    obs_im = double(obs_im);
  end

  if (INTENSITY_SCALING)
    obs_im = obs_im.*INTENSITY_SCALING;
  end

  if (EXTRA_BLUR)
    f = fspecial('gaussian',[1 1]*EXTRA_BLUR_SIZE,EXTRA_BLUR_SIGMA);
    for c=1:obs_imz
      obs_im(:,:,c) = conv2(obs_im(:,:,c),f,'same');
    end
  end
  
  if (SATURATION_MASK)
   
    %%% Use intensity channel
    sat = (obs_im(:,:,1) > SATURATION_THRESHOLD);
    q = conv2(double(sat),ones(size(blur_kernel)),'same');
    mask = (q>0);
  
  else
    mask = zeros(size(obs_im,1),size(obs_im,2));
  end


  %% Images are linear so 
  if AUTOMATIC_PATCH
    [obs_im_tmp,PATCH_LOCATION] = automatic_patch_selector(obs_im.^(1/DEFAULT_GAMMA),max(PATCH_SIZE),AUTOMATIC_PATCH_CENTER_WEIGHT,mask);
  end
    
  if strcmp(CAMERA_TYPE,'unknown')
    %%% leave [0:255] alone
  else
    %%% load up map for file
    load(CAMERA_TYPE);
    im_c = obs_im(:);
  
    if COLOR
      
     for a=1:3 %% each color plane
        new_col = interp1q([0:255]',response_curves(a,:)',im_c);
        obs_im(:,:,a) = reshape(new_col,obs_imy,obs_imx);
      end
      
    else
      
      new_col = interp1q([0:255]',response_curves(4,:)',im_c);
      obs_im = reshape(new_col,obs_imy,obs_imx);
   
    end
    
  end

  if COLOR %%% Convert to YIQ colorspace
    obs_im = rgb2ntsc(obs_im);
  end  
       
  
  % final level is high res version
  obs_im_all = obs_im;
  mask_all   = mask;
  
  if ~strcmp(GRADIENT_MODE,'steer')
    
    % select gradient filter
    if strcmp(GRADIENT_MODE,'haar')
      kx = [1 -1]; ky = [1 -1]';
    elseif strcmp(GRADIENT_MODE,'laplace')
      kx = [1 -2 1]; ky = [1 -2 1]';
    else 
      error foo
    end
  
    for c=1:obs_imz
      obs_grad_all_x(:,:,c) = conv2(obs_im_all(:,:,c),kx,'valid');
      obs_grad_all_y(:,:,c) = conv2(obs_im_all(:,:,c),ky,'valid'); 
      
      %% just take x for the time being
      yy = min(size(obs_grad_all_x,1),size(obs_grad_all_y,1));
      xx = min(size(obs_grad_all_x,2),size(obs_grad_all_y,2));
      
      obs_grad_all(:,:,c) = [obs_grad_all_x(1:yy,1:xx,c),obs_grad_all_y(1:yy,1:xx,c)];
    end
    
  end
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chop out patch and form scale pyramid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % Chop out patch in intensity and gradient space
  py = 0; 
  px = 0;
  obs_im = obs_im_all(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px,:);

  mask = mask_all(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px);
    
  if ((any(findstr(RESIZE_MODE,'matlab')) & ~strcmp(GRADIENT_MODE,'steer')) | (RESIZE_STEP~=2))
    
    %% chop out gradient patch too
    obs_grad_x = obs_grad_all_x(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px,:);
    obs_grad_y = obs_grad_all_y(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px,:);    

       
    if RESCALE_THEN_GRAD
      %%%% rescale images then take gradient 
      
      % use matlab's imresize function
      for s = 2:NUM_SCALES
        if strcmp(RESIZE_MODE,'matlab_nearest')
          obs_im_scale{NUM_SCALES-s+1} = imresize(obs_im,(1/RESIZE_STEP)^(s-1),'nearest');
        elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
          obs_im_scale{NUM_SCALES-s+1} = imresize(obs_im,(1/RESIZE_STEP)^(s-1),'bilinear');       
        elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
          obs_im_scale{NUM_SCALES-s+1} = imresize(obs_im,(1/RESIZE_STEP)^(s-1),'bicubic');       
        else
          error('Must use Matlab imresize funcition if RESIZE_STEP~=2');
        end
   
        %%% Compute saturation mask
        mask_scale{NUM_SCALES-s+1} = ceil(abs(imresize(mask,(1/RESIZE_STEP)^(s-1),'nearest')));
        mask_scale{NUM_SCALES-s+1} = mask_scale{NUM_SCALES-s+1}(2:end-1,2:end-1);
        
      end

      
      % final level is high res version
      obs_im_scale{NUM_SCALES} = obs_im;
        mask_scale{NUM_SCALES} = mask(2:end-1,2:end-1);
    
      % select gradient filter
      if strcmp(GRADIENT_MODE,'haar')
        kx = [0 1 -1]; ky = [0 1 -1]';
      elseif strcmp(GRADIENT_MODE,'laplace')
        kx = [1 -2 1]; ky = [1 -2 1]';
      else 
        error foo
      end
      
      % apply to resized images
      for s = 1:NUM_SCALES
        for c = 1:obs_imz
          obs_grad_scale_x{s}(:,:,c) = conv2(obs_im_scale{s}(:,:,c),kx,'valid');
          obs_grad_scale_y{s}(:,:,c) = conv2(obs_im_scale{s}(:,:,c),ky,'valid'); 
        end   
      end
      
    else %%% GRAD THEN RESCALE

      % final level is high res version
      obs_grad_scale_x{NUM_SCALES} = obs_grad_x;
      obs_grad_scale_y{NUM_SCALES} = obs_grad_y;
      mask_scale{NUM_SCALES}       = mask;
       
      if MANUAL_KERNEL
        manual_kernel{NUM_SCALES} = approximate_kernel; %%% loaded up in script
      end
      
      % use matlab's imresize function
      for s = 2:NUM_SCALES
        if strcmp(RESIZE_MODE,'matlab_nearest')
          obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_x,(1/RESIZE_STEP)^(s-1),'nearest');
          obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_y,(1/RESIZE_STEP)^(s-1),'nearest');
        elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
          obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_x,(1/RESIZE_STEP)^(s-1),'bilinear');       
          obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_y,(1/RESIZE_STEP)^(s-1),'bilinear');       
         %obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_scale_x{NUM_SCALES-s+2},(1/RESIZE_STEP),'bilinear');       
         %obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_scale_y{NUM_SCALES-s+2},(1/RESIZE_STEP),'bilinear');       
              
        elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
          obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_x,(1/RESIZE_STEP)^(s-1),'bicubic');       
          obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_y,(1/RESIZE_STEP)^(s-1),'bicubic');       
        else
          error('Must use Matlab imresize funcition if RESIZE_STEP~=2');
        end

        %%% Compute saturation mask
        mask_scale{NUM_SCALES-s+1} = ceil(abs(imresize(mask,(1/RESIZE_STEP)^(s-1),'nearest')));
       
        %%% Add in additional blur for all level except final
        %obs_grad_scale_x{NUM_SCALES-s+1} = conv2(obs_grad_scale_x{NUM_SCALES-s+1},fspecial('gaussian',[5 5],0.5),'same');
        %%% Add in additional blur for all level except final
        %obs_grad_scale_y{NUM_SCALES-s+1} = conv2(obs_grad_scale_y{NUM_SCALES-s+1},fspecial('gaussian',[5 5],0.5),'same');
        %fprintf('Adding in 0.5 blur to each level\n');
      end
     
   
      
    end
    
  else %%% NOT USING MATLAB
    %%% NOT SURE IF IT WILL WORK IN COLOR MODE
    %%% these do both rescaling and gradients all together...
    %%% use wavelet/pyramid code to do scaling
 
    if strcmp(GRADIENT_MODE,'haar')
      [pyr,ind] = buildWpyr(double(obs_im),NUM_SCALES,'haar');
      for s = 1:NUM_SCALES
        obs_grad_scale_x{s} = wpyrBand(pyr,ind,NUM_SCALES-s+1,1);
        obs_grad_scale_y{s} = wpyrBand(pyr,ind,NUM_SCALES-s+1,2);  
      end 
    elseif strcmp(GRADIENT_MODE,'steer')
      [pyr,ind] = buildSpyr(double(obs_im),NUM_SCALES);
      for s = 1:NUM_SCALES
        obs_grad_scale_x{s} = spyrBand(pyr,ind,NUM_SCALES-s+1,1);
        obs_grad_scale_y{s} = spyrBand(pyr,ind,NUM_SCALES-s+1,2);  
      end   
    else
      error foo
    end
  
  end
  
  % Just take x-gradients only at the moment
  %obs_grad_scale = [obs_grad_scale_x]; 
  
  
  if (RESCALE_THEN_GRAD)
  
    for s = 1:NUM_SCALES
      obs_grad_scale{s} = [ obs_grad_scale_x{s}(2:end-1,:) , obs_grad_scale_y{s}(:,2:end-1) ];
    end
  
    
  else
    
    %%% Concatenate the gradient images together
    for s = 1:NUM_SCALES
      obs_grad_scale{s} = [ obs_grad_scale_x{s} , obs_grad_scale_y{s} ];
    end
  
  end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess blur kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(findstr(RESIZE_MODE,'matlab')) | (RESIZE_STEP~=2)
 
  blur_kernel_scale{NUM_SCALES} = blur_kernel;
  
  %%% use matlab's imresize function
  for s = 2:NUM_SCALES

    dims = size(blur_kernel_scale{NUM_SCALES}) * (1/RESIZE_STEP)^(s-1);
    dims = dims + (1-mod(dims,2)); %% make odd size

  
    if (min(dims)<4)
      
      %% Blur kernel first, then resize (since it will use nearest neighbour)
      h = fspecial('gaussian',dims,1);
      blur_kernel_scale{NUM_SCALES-s+1} = imresize(conv2(blur_kernel,h),dims,'nearest');
   
    else
      
      if strcmp(RESIZE_MODE,'matlab_nearest')
        blur_kernel_scale{NUM_SCALES-s+1} = imresize(blur_kernel,dims,'nearest');
      elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
        blur_kernel_scale{NUM_SCALES-s+1} = imresize(blur_kernel,dims,'bilinear');       
      elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
        blur_kernel_scale{NUM_SCALES-s+1} = imresize(blur_kernel,dims,'bicubic');       
      else
        error foo
      end
    
    end
    
    blur_kernel_scale{NUM_SCALES-s+1} = blur_kernel_scale{NUM_SCALES-s+1} / sum(blur_kernel_scale{NUM_SCALES-s+1}(:));

  end
  
    
else
  
  %%% use wavelet/pyramid code to do scaling
  
  blur_kernel_scale{NUM_SCALES} = blur_kernel;
  for s = 2:NUM_SCALES
    blur_kernel_scale{NUM_SCALES-s+1} = blurDn(blur_kernel_scale{NUM_SCALES-s+2},1,RESIZE_MODE);
    blur_kernel_scale{NUM_SCALES-s+1} = blur_kernel_scale{NUM_SCALES-s+1} / sum(blur_kernel_scale{NUM_SCALES-s+1}(:));
  end 
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate blurred image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SYNTHETIC
  
  % the obs image is really the sharp one if we are synthetic
  % so copy it to true_im
  true_im_all = obs_im_all;
  true_grad_all_x = obs_grad_all_x;
  true_grad_all_y = obs_grad_all_y;
  true_grad_all = obs_grad_all;
  true_im = obs_im; % the patch
  true_grad_scale = obs_grad_scale; % gradient pyramid
  true_grad_scale_x = obs_grad_scale_x; % x-grad pyr
  true_grad_scale_y = obs_grad_scale_y; % y-grad pyr
  
  % now blur the obs_im

  if (FFT_MODE)
    obs_im_all = real(ifft2(fft2(obs_im_all) .* fft2(blur_kernel,size(obs_im_all,1),size(obs_im_all,2))));
  else
    obs_im_all = conv2(obs_im_all,blur_kernel,'same');
  end
  
  %%% Now round back to unit8 image
  %obs_im_all = double(uint8(obs_im_all));
   
  
if COLOR
  %% TODO
else

    % Chop out patch in intensity and gradient space
    py = 0; %round((floor(length(ky/2)) + floor(size(blur_kernel,1)))/2);
    px = 0; %round((floor(length(kx/2)) + floor(size(blur_kernel,2)))/2);
    obs_im = obs_im_all(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px,:);
     
   
    if ((any(findstr(RESIZE_MODE,'matlab')) | ~strcmp(GRADIENT_MODE,'steer'))   | (RESIZE_STEP~=2))
 
    if RESCALE_THEN_GRAD
       %%%% rescale images then take gradient 
      
       % use matlab's imresize function
      for s = 2:NUM_SCALES
        if strcmp(RESIZE_MODE,'matlab_nearest')
          obs_im_scale{NUM_SCALES-s+1} = imresize(obs_im,(1/RESIZE_STEP)^(s-1),'nearest');
        elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
          obs_im_scale{NUM_SCALES-s+1} = imresize(obs_im,(1/RESIZE_STEP)^(s-1),'bilinear');       
        elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
          obs_im_scale{NUM_SCALES-s+1} = imresize(obs_im,(1/RESIZE_STEP)^(s-1),'bicubic');       
        else
          error('Must use Matlab imresize funcition if RESIZE_STEP~=2');
        end
      end
      
      % final level is high res version
      obs_im_scale{NUM_SCALES} = obs_im;
      
      % select gradient filter
      if strcmp(GRADIENT_MODE,'haar')
        kx = [1 -1]; ky = [1 -1]';
      elseif strcmp(GRADIENT_MODE,'laplace')
        kx = [1 -2 1]; ky = [1 -2 1]';
      else 
        error foo
      end
      
      % apply to resized images
      for s = 1:NUM_SCALES
        obs_grad_scale_x{s} = conv2(obs_im_scale{s},kx,'valid');
        obs_grad_scale_y{s} = conv2(obs_im_scale{s},ky,'valid'); 
      end   
      
    else %%% grad then rescale
    
  
      %% take gradients
      if ~strcmp(GRADIENT_MODE,'steer')
        % select gradient filter
        if strcmp(GRADIENT_MODE,'haar')
          kx = [1 -1]; ky = [1 -1]';
        elseif strcmp(GRADIENT_MODE,'laplace')
          kx = [1 -2 1]; ky = [1 -2 1]';
        else 
          error foo
        end
        
        %% take gradients of whole image
        obs_grad_all_x = conv2(obs_im_all,kx,'valid');
        obs_grad_all_y = conv2(obs_im_all,ky,'valid'); 
 
        %% just take x for the time being
        yy = min(size(obs_grad_all_x,1),size(obs_grad_all_y,1));
        xx = min(size(obs_grad_all_x,2),size(obs_grad_all_y,2));
    
        obs_grad_all = [obs_grad_all_x(1:yy,1:xx),obs_grad_all_y(1:yy,1:xx)];
     
      end

      % chop out patch from gradients of burred image
      obs_grad_x = obs_grad_all_x(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px,:);
      obs_grad_y = obs_grad_all_y(PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py,PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px,:);    

      % use matlab's imresize function
      for s = 2:NUM_SCALES
        if strcmp(RESIZE_MODE,'matlab_nearest')
          obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_x,(1/RESIZE_STEP)^(s-1),'nearest');
          obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_y,(1/RESIZE_STEP)^(s-1),'nearest');
        elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
          obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_x,(1/RESIZE_STEP)^(s-1),'bilinear');       
          obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_y,(1/RESIZE_STEP)^(s-1),'bilinear');       
        elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
          obs_grad_scale_x{NUM_SCALES-s+1} = imresize(obs_grad_x,(1/RESIZE_STEP)^(s-1),'bicubic');       
          obs_grad_scale_y{NUM_SCALES-s+1} = imresize(obs_grad_y,(1/RESIZE_STEP)^(s-1),'bicubic');       
        else
          error('Must use Matlab imresize funcition if RESIZE_STEP~=2');
        end
      end
      
      % final level is high res version
      obs_grad_scale_x{NUM_SCALES} = obs_grad_x;
      obs_grad_scale_y{NUM_SCALES} = obs_grad_y;
      
    end
    
  else 
    
    %%% use wavelet/pyramid code to do scaling
    
    if strcmp(GRADIENT_MODE,'haar')
      [pyr,ind] = buildWpyr(double(obs_im),NUM_SCALES,'haar');
      for s = 1:NUM_SCALES
        obs_grad_scale_x{s} = wpyrBand(pyr,ind,NUM_SCALES-s+1,1);
        obs_grad_scale_y{s} = wpyrBand(pyr,ind,NUM_SCALES-s+1,2);  
      end 
    elseif strcmp(GRADIENT_MODE,'steer')
      [pyr,ind] = buildSpyr(double(obs_im),NUM_SCALES);
      for s = 1:NUM_SCALES
        obs_grad_scale_x{s} = spyrBand(pyr,ind,NUM_SCALES-s+1,1);
        obs_grad_scale_y{s} = spyrBand(pyr,ind,NUM_SCALES-s+1,2);  
      end   
    else
      error foo
    end
  
  end
  
  % Just take x-gradients only at the moment
  %obs_grad_scale = [obs_grad_scale_x]; 
  
  %%% Concatenate the gradient images together
  for s = 1:NUM_SCALES
    obs_grad_scale{s} = [ obs_grad_scale_x{s} , obs_grad_scale_y{s} ];
  end
  
end

else %%% NOT SYNTHIETC

  if (FFT_MODE)
    for s=1:NUM_SCALES
      
      obs_grad_scale_old{s} = obs_grad_scale{s};
      db = delta_kernel(size(blur_kernel_scale{s},1));
      for c=1:obs_imz
        obs_grad_scale{s}(:,:,c) = real(ifft2(fft2(obs_grad_scale{s}(:,:,c)).*fft2(db,size(obs_grad_scale{s},1),size(obs_grad_scale{s},2))));
      
      end
 
      %%% translate mask too...
      mask_scale{s} = conv2(double(mask_scale{s}),db,'same');
 
    end
  end

end  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get sizes of everything and make masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1:NUM_SCALES

  % blur 
  K(s) = size(blur_kernel_scale{s},1);
  L(s) = size(blur_kernel_scale{s},2);
  hK(s) = floor(K(s)/2); hL(s) = floor(L(s)/2);

  % image 
  M(s) = size(obs_grad_scale{s},1);
  N(s) = size(obs_grad_scale{s},2);

  if BLUR_MASK
    for a=1:BLUR_COMPONENTS
      [xx,yy] = meshgrid(-hK(s):hK(s),-hL(s):hL(s));
      spatial_blur_mask{s}(a,:) = log(normMDpdf([xx(:)';yy(:)'],[0 0]',eye(2)*K(s)*BLUR_MASK_VARIANCES(a)));
    end
  else
    spatial_blur_mask{s} = zeros(BLUR_COMPONENTS,K(s)*L(s));
  end
 
  if (SATURATION_MASK==2)
    spatial_image_mask{s} = [mask_scale{s} * SATURATION_PRESCISION,mask_scale{s} * SATURATION_PRESCISION];
  else
    spatial_image_mask{s} = [zeros(size(mask_scale{s})),zeros(size(mask_scale{s}))];
  end
  
  % observed
  if (FFT_MODE)
    I(s) = 2*M(s);
    J(s) = 2*N(s);
    Dp{s} = zeros(I(s),J(s));
  
    Dp{s}(K(s):M(s),L(s):N(s)/2) = 1; %% x-plane
    Dp{s}(K(s):M(s),L(s)+N(s)/2:N(s)) = 1; %% y-plane
    
    D{s} = padarray(obs_grad_scale{s},[M(s) N(s)],0,'post');

    %%% Now add in saturation mask
    if (SATURATION_MASK==1)
      Dp{s} = Dp{s} .* padarray(1-[mask_scale{s},mask_scale{s}],[M(s) N(s)],0,'post');
    end
 
  else
    I(s) = M(s);
    J(s) = N(s);
   
    D{s} = obs_grad_scale{s}; 
     % make masks
   
    for c=1:obs_imz
      Dp{s}(:,:,c) = zeros(I(s),J(s));
     
      Dp{s}(hK(s)+1:M(s)-hK(s),hL(s)+1:N(s)/2-hL(s),c) = 1;
      Dp{s}(hK(s)+1:M(s)-hK(s),N(s)/2+hL(s)+1:N(s)-hL(s),c) = 1;
      
      %%% Now add in saturation mask
      if (SATURATION_MASK==1)
        Dp{s}(:,:,c) = Dp{s}(:,:,c) .* (1-[mask_scale{s},mask_scale{s}]);
      end
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:NUM_SCALES
  
   
  if ((NUM_SCALES-s+1)>8)
    priors(NUM_SCALES-s+1) = priors(8);
  end    
  
  % Loop over scales
  % nDims , size_per_dims, #compoents in prior, #prior type,lock prior,update]            
  dimensions = [1       1         1                0           0         1;
                1       K(s)*L(s) BLUR_COMPONENTS  BLUR_PRIOR  BLUR_LOCK 1;
                obs_imz M(s)*N(s) IMAGE_COMPONENTS IMAGE_PRIOR 1         1];   

 % Initialize
  if (s==1) % first iteration
  
    if SYNTHETIC
      [x1,x2] = initialize_parameters2(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s},true_grad_scale{s},INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,FIRST_INIT_MODE_IMAGE,FIRST_INIT_MODE_BLUR,obs_im,blur_kernel_scale{NUM_SCALES},spatial_image_mask{s},priors(NUM_SCALES-s+1),FFT_MODE,COLOR,1);
    else
      [x1,x2] = initialize_parameters2(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,FIRST_INIT_MODE_IMAGE,FIRST_INIT_MODE_BLUR,obs_im,blur_kernel_scale{NUM_SCALES},spatial_image_mask{s},priors(NUM_SCALES-s+1),FFT_MODE,COLOR,1);
   end
%  elseif (s==2) % 2nd iterations
%    if SYNTHETIC
%      [x1,x2] = initialize_parameters2(obs_grad_scale{s},[],[],blur_kernel_scale{s},true_grad_scale{s},INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,FIRST_INIT_MODE_IMAGE,FIRST_INIT_MODE_BLUR);
%    else
%      [x1,x2] = initialize_parameters2(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,INIT_MODE_IMAGE,'true');
%   end
  else % as we work up pyramid

    if SYNTHETIC
      [x1,x2] = initialize_parameters2(obs_grad_scale{s},new_blur{s-1},new_grad{s-1},blur_kernel_scale{s},true_grad_scale{s},INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,INIT_MODE_IMAGE,INIT_MODE_BLUR,obs_im,blur_kernel_scale{NUM_SCALES},spatial_image_mask{s},priors(NUM_SCALES-s+1),FFT_MODE,COLOR,1);
    else
      [x1,x2] = initialize_parameters2(obs_grad_scale{s},new_blur{s-1},new_grad{s-1},[],[],INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,INIT_MODE_IMAGE,INIT_MODE_BLUR,obs_im,blur_kernel_scale{NUM_SCALES},spatial_image_mask{s},priors(NUM_SCALES-s+1),FFT_MODE,COLOR,1);     
    end
  end

  m = x1./x2;
  mx_init{s} = reshape(m(2+K(s)*L(s):end),M(s),N(s),obs_imz);
  me_init{s} = reshape(m(2:1+K(s)*L(s)),K(s),L(s));
 
  % Call routine
  [ensemble,D_log{s},gamma_log]=train_ensemble_main6(dimensions,x1,x2,'',['Scale=' int2str(s)],[CONVERGENCE 0 NOISE_INIT 0 0 MAX_ITERATIONS 0],D{s},Dp{s},I(s),J(s),K(s),L(s),M(s),N(s),priors(NUM_SCALES-s+1),FFT_MODE,spatial_blur_mask{s},1-(spatial_image_mask{s}(:)>0)); 
 
  
  % Extract mean blur/image
  me_est{s} =reshape(train_ensemble_get(2,dimensions,ensemble.mx),K(s),L(s));
  me_est2{s} =reshape(train_ensemble_get(2,dimensions,ensemble.mx2),K(s),L(s));

  mx_est{s} =reshape(train_ensemble_get(3,dimensions,ensemble.mx),M(s),N(s),obs_imz);
  mx_est2{s} =reshape(train_ensemble_get(3,dimensions,ensemble.mx2),M(s),N(s),obs_imz);
  
  
  if (s~=NUM_SCALES)

    % Use solution for next level up...
    [new_grad{s},new_blur{s}] = move_level(mx_est{s},me_est{s},K(s+1),L(s+1),M(s+1),N(s+1),UPSAMPLE_MODE,RESIZE_STEP,CENTER_BLUR);
  
    if MANUAL_KERNEL
      % Overwrite upsampled version with manual version
      new_blur{s} = blur_kernel_scale{s+1};
    end
        
  end

  %%% save everything out

  save(['tmp_',char(CONFIG_FNAME)],'mx_init','me_init','mx_est','me_est','new_grad','new_blur','ensemble','D_log');
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct image patch to intensity space from gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% First try reconstructing image from gradients (Yair's code)
%%% Note that we can only do this for the image patch - not the whole image
if COLOR
  for c=1:obs_imz
    final_dx = mx_est{NUM_SCALES}(:,1:N(NUM_SCALES)/2,c);
    final_dy = mx_est{NUM_SCALES}(:,N(NUM_SCALES)/2+1:N(NUM_SCALES),c);
    final_dx(1,:) = 0; final_dx(:,1) = 0;
    final_dx(end,:) = 0; final_dx(:,end) = 0;
    final_dy(1,:) = 0; final_dy(:,1) = 0;
    final_dy(end,:) = 0; final_dy(:,end) = 0;
    [obs_im_recon(:,:,c),k] = reconsEdge3(final_dx,final_dy); 
  end
else
  final_dx = mx_est{NUM_SCALES}(:,1:N(NUM_SCALES)/2);
  final_dy = mx_est{NUM_SCALES}(:,N(NUM_SCALES)/2+1:N(NUM_SCALES));
  final_dx(1,:) = 0; final_dx(:,1) = 0;
  final_dx(end,:) = 0; final_dx(:,end) = 0;
  final_dy(1,:) = 0; final_dy(:,1) = 0;
  final_dy(end,:) = 0; final_dy(:,end) = 0;
  [obs_im_recon,k] = reconsEdge3(final_dx,final_dy); 
end


CONSISTENCY_GRAD_ITS = 0;
for a=1:CONSISTENCY_GRAD_ITS

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Rerun infernce using estimated intensities.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% recompute gradients
  kx = [1 -1]; ky = [1 -1]';
  
  for c=1:obs_imz
    obs_grad_recon_x(:,:,c) = conv2(obs_im_recon(:,:,c),kx,'same');
    obs_grad_recon_y(:,:,c) = conv2(obs_im_recon(:,:,c),ky,'same'); 
       
    obs_grad_recon(:,:,c) = [obs_grad_recon_x(:,:,c),obs_grad_recon_y(:,:,c)];
  end
  
  %% Run infernce again
  [x1,x2] = initialize_parameters2(obs_grad_scale{end},me_est{end},obs_grad_recon,[],[],INIT_PRESCISION,IMAGE_PRIOR,IMAGE_COMPONENTS,INIT_MODE_IMAGE,INIT_MODE_BLUR,obs_im,blur_kernel_scale{NUM_SCALES},spatial_image_mask{end},priors(1),FFT_MODE,COLOR,1);   
  
  [ensemble,D_log_recon,gamma_log]=train_ensemble_main6(dimensions,x1,x2,'',['Final scale refinement'],[CONVERGENCE 0 NOISE_INIT 0 0 MAX_ITERATIONS 0],D{end},Dp{end},I(end),J(end),K(end),L(end),M(end),N(end),priors(1),FFT_MODE,spatial_blur_mask{end},1-(spatial_image_mask{end}(:)>0)); 

  me_est{end} =reshape(train_ensemble_get(2,dimensions,ensemble.mx),K(end),L(end));
  mx_est{end} =reshape(train_ensemble_get(3,dimensions,ensemble.mx),M(end),N(end),obs_imz);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Reconstruct image patch again
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% First try reconstructing image from gradients (Yair's code)
  %%% Note that we can only do this for the image patch - not the whole image
  if COLOR
    for c=1:obs_imz
      final_dx = mx_est{NUM_SCALES}(:,1:N(NUM_SCALES)/2,c);
      final_dy = mx_est{NUM_SCALES}(:,N(NUM_SCALES)/2+1:N(NUM_SCALES),c);
      final_dx(1,:) = 0; final_dx(:,1) = 0;
      final_dx(end,:) = 0; final_dx(:,end) = 0;
      final_dy(1,:) = 0; final_dy(:,1) = 0;
      final_dy(end,:) = 0; final_dy(:,end) = 0;
      [obs_im_recon(:,:,c),k] = reconsEdge3(final_dx,final_dy); 
    end
  else
    final_dx = mx_est{NUM_SCALES}(:,1:N(NUM_SCALES)/2);
    final_dy = mx_est{NUM_SCALES}(:,N(NUM_SCALES)/2+1:N(NUM_SCALES));
    final_dx(1,:) = 0; final_dx(:,1) = 0;
    final_dx(end,:) = 0; final_dx(:,end) = 0;
    final_dy(1,:) = 0; final_dy(:,1) = 0;
    final_dy(end,:) = 0; final_dy(:,end) = 0;
    [obs_im_recon,k] = reconsEdge3(final_dx,final_dy); 
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save to mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(char(CONFIG_FNAME),'ensemble','dimensions','obs_im_recon','obs_im_orig','obs_im_all','obs_grad_scale','mask_all','blur_kernel_scale','obs_im_recon','blur_kernel','mx_init','me_init','mx_est','me_est','new_grad','new_blur','PATCH_LOCATION','-append');
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now do actual deblurring of whole image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(IMAGE_RECONSTRUCTION,'lucy')

  %%% Run RL script
  [deblurred_im,kernel_out] = fiddle_lucy3(char(CONFIG_FNAME),LUCY_ITS,1,SCALE_OFFSET,KERNEL_THRESHOLD);
  
else
  error('Currently only implemented for Richardson-Lucy');
end

% save final image and kernel
save(char(CONFIG_FNAME),'deblurred_im','kernel_out','-append');
