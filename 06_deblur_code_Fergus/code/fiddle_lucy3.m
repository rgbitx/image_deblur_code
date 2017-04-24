function [out,blur_kernel]=fiddle_lucy3(file_name,LUCY_ITS_IN,SAVE_TO_DISK,SCALE_OFFSET_IN,THRESHOLD_IN)
%
% Rountine to call Richardson-Lucy algorithm after kernel inference  
%
% Inputs:
% LUCY_ITS - integer. Default = 10. Sensible values 5 to
%    30. Number of Lucy-Richardson iterations to use when unblurring
%    the image using the inferred kernel. 10 is the default but if the
%    kernel is very long and thin you might need more. If you turn it
%    up too much it amplifies the noise in the image.
%
% SAVE_TO_DISK - binary. Turns saving of deblurred image and final kernel
% (post thresholding) to disk or not.
%
%
% SCALE_OFFSET - integer. Default = 0. Sometimes it may not be
%      possible to deblur the image at full resolution due to high
%      noise levels (see discussion in 7.3). If this is the case, you
%      can tell the fiddle_lucy3 function to use a coarser scale. The
%      value of SCALE_OFFSET dicates how many scale levels you drop
%      down before selecting the kernel. i.e. SCALE_OFFSET = 1 will use
%      a kernel that is one scale level (RESIZE_STEP smaller) than the
%      full resolution kernel.
%
% THRESHOLD - float. Default = 7, sensible range is 5 to
%    15. This is a threshold on kernel intensities, applied after the
%    whole inference stage, designed to remove noise from the
%    kernel. It is a dynamic thresold which is the percentage of the
%    maximum value of the kernel. Since it is a bit dependent on the
%    intensity profile of the kernel, some manual tuning may be needed
%    to get the optimal value. If you don't want to use the threshold at all, set it to a
%    large number (i.e. 100). If you set it too high, it will start to
%    erode the structure of the kernel. This parameter is a bit of a
%    hack - in reality you shouldn't need it, but it can make quite a
%    bit of difference if the inferred kernel is noisy. 
%
% Outputs:
% 1. out - Deblurred image  
% 2. blur_kernel - blur kernel after thresholding step.

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology
  
SHOW_BOX = 1; % show gray rectangle on blurry image indicating selected region
EDGE_CROP = 1; % crop edge of deblurred image off (since you can't
               % recover it).

if (nargin==4)
  THRESHOLD_IN = 10;
end

if (nargin==3)
  SCALE_OFFSET_IN = 0;
  THRESHOLD_IN = 10;
end

if (nargin==2)
  SCALE_OFFSET_IN = 0;
  SAVE_TO_DISK_IN = 0;
  THRESHOLD_IN = 10;
end

if (nargin==1)
  LUCY_ITS_IN = 10;
  SCALE_OFFSET_IN = 0;
  SAVE_TO_DISK_IN = 0;
  THRESHOLD_IN = 10;
end  

%%% Load up model
load(file_name);

%% override values in save file
SCALE_OFFSET = SCALE_OFFSET_IN;
LUCY_ITS = LUCY_ITS_IN;
THRESHOLD = THRESHOLD_IN;

%%% Get blur kernel
blur_kernel = me_est{end-SCALE_OFFSET}/sum(me_est{end-SCALE_OFFSET}(:));

%%% Threshold kernel
threshold = max(blur_kernel(:))/THRESHOLD;
z=find(blur_kernel(:)<threshold);
blur_kernel(z)=0;

% check kernel sums to 1
blur_kernel = blur_kernel / sum(blur_kernel(:));

%%% Resize image
if (PRESCALE~=1)
  obs_im = imresize(obs_im,PRESCALE,'bilinear');
end

%%% Rescale if not using final scale
if (SCALE_OFFSET>0)
  obs_im = imresize(obs_im,(1/sqrt(2))^SCALE_OFFSET,'bilinear');
end

%%% Take out gamma corection
if (GAMMA_CORRECTION~=1)
  obs_im_gam = ((double(obs_im).^GAMMA_CORRECTION)/(256^(GAMMA_CORRECTION-1)));
else
  obs_im_gam = obs_im;
end

fprintf('.');

% use edgetaper to sort out edge effects
% obs_im_gam = edgetaper(obs_im_gam,blur_kernel); 

fprintf('.');

% run RL
out = deconvlucy(obs_im_gam,blur_kernel,LUCY_ITS);  

fprintf('.');

     
if (GAMMA_CORRECTION~=1)
  %%% Put gamma back in
  out = double(out).^(1/GAMMA_CORRECTION);
else
  out = double(out);
end

%% shift and rescale to make 0 to 1 image
out = out - min(out(:));
out = out / max(out(:));

%% now do histogram equalization
out = histmatch(out,uint8(obs_im));

%% apply box to im_col
if SHOW_BOX
  aa = round(AXIS * (1/sqrt(2))^SCALE_OFFSET);
  obs_im(aa(3),aa(1):aa(2),:) = 75;
  obs_im(aa(4),aa(1):aa(2),:) = 75;
  obs_im(aa(3):aa(4),aa(1),:) = 75;
  obs_im(aa(3):aa(4),aa(2),:) = 75;
end

%%% Crop image to avoid edge artifacts
if EDGE_CROP
  edge_offset = floor(size(blur_kernel,1)/2);
  out = out(edge_offset+1:end-edge_offset-1,edge_offset+1:end-edge_offset-1,:);
  obs_im = obs_im(edge_offset+1:end-edge_offset-1,edge_offset+1:end-edge_offset-1,:);
end

%%%%% Show images/kernel

%% Plot output image
figure; imagesc(out); title('Output'); axis equal;
%% Plot original, blurry image
figure; imagesc(obs_im); title('Original image'); axis equal;

%%% Plot kernels
h=figure; imagesc(blur_kernel); colormap(gray);
axis square; title('Inferred kernel');

if SAVE_TO_DISK
  
  % save kernel
  ExportFig(h,[file_name,'_kernel']);
  % save blurry
  imwrite(uint8(obs_im),[file_name,'_blurry.jpg'],'jpg','Quality',100);
  % save blurry
  imwrite(uint8(out),[file_name,'_final.jpg'],'jpg','Quality',100);

end
