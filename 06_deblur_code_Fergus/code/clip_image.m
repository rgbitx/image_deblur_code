function [im] = clip_image(im,minval,maxval)

% Author: Bryan Russell
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% [im] = clip_image(im,minval,maxval)
% 
% Description:
% Clips image to have values between minval and maxval (i.e. minval <= im
% <= maxval).  
%
% Inputs:
% im - The image.
% minval - Minimum value.
% maxval - Maximum value.
%
% Outputs:
% im - Clipped image
%
  im(im < minval) = minval;
  im(im > maxval) = maxval;
