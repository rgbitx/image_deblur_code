function [out_im,patch_location]=automatic_patch_selector(im,patch_size,weight,sat_mask)
  
% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

  SMOOTH_SIGMA = 3;
  
  %% Assume input image is 0:255
  %% patch_size is integer (odd)
  %% sat_mask is binary & same size as im
    
  %% weight is tuning parameter btw. variance and non-saturated pixels
    
  %%% Get size of input im
  [II,JJ] = size(im);
  
  %%% Compute centre weighting mask
  [xx,yy] = meshgrid([1:JJ]-round(JJ/2),[1:II]-round(II/2));
  centre_weight_mask = exp(-weight/(JJ^2)*(xx.^2+yy.^2));
  II = II*2; JJ = JJ*2;
  %% shift by patch_size
  centre_weight_mask = real( ifft2( fft2(centre_weight_mask,II,JJ) .* fft2(delta_kernel(patch_size),II,JJ) ) );
  
  %%% Get patch mask
  pmask = ones(patch_size)/patch_size.^2;
  
  %%% Find patch with largest variance 
  ei2 = real( ifft2( fft2(im.^2,II,JJ) .* fft2(pmask,II,JJ) ) );
  mu2 = real( ifft2( fft2(im,II,JJ) .* fft2(pmask,II,JJ) ) ).^2;
  w = ei2 - mu2;
  
  %%% Compute convolution with sautration mask
  q = real( ifft2( fft2((sat_mask),II,JJ) .* fft2(pmask,II,JJ) ) );
  %% q is small if more pixels are available for use 


  combined = (centre_weight_mask).*w./(q*mean(im(:)).^2+1); %% more variance, less saturation

  %% now find stable maximum (smooth resonse image)
  f = fspecial('gaussian',[8 8],SMOOTH_SIGMA);
  combined = real( ifft2( fft2(combined,II,JJ) .* fft2(f,II,JJ) ) );

  %%% crop to avoid edge effects
  combined = combined(patch_size:II/2,patch_size:JJ/2);
  
  %%% find max
  [tmp,mm] = max(combined(:));
  [sy,sx] = ind2sub(size(combined),mm);

  %%% get coords for axis
  patch_location = [sx sy] -1;
  
  %% chop out patch
  out_im = im(sy-1:sy-2+patch_size,sx-1:sx-2+patch_size);
  
  
  
