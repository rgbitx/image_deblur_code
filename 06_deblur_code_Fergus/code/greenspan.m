function [en,L0] = greenspan(im,S)

% Author: Bryan Russell
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% GREENSPAN - Uses the algorithm in Greenspan,Anderson,Akber, "Image
% Enhancement by Nonlinear Extrapolation in Frequency Space" for images.
%
% [EN,L0] = GREENSPAN(IM,S) enhances the image IM using the settings
% specified in S and returns the enhanced, super-resolved output image EN
% and the inferred high frequency information L0.
%
% See also CREATE_GREENSPAN_SETTINGS
%
% <brussell@csail.mit.edu>
  
  z = 2^(S.factor);
  L1 = im-rconv2(im,S.lo_filt);
  L0 = upConv(L1,z^2*S.lo_filt,'reflect1',[z z],[1 1],z*size(im));
% $$$   L0 = upConv(L1,4*S.lo_filt,'reflect1',[2 2],[1 1],2*size(im));
  L0 = S.s*clip_image(L0,-(1-S.c)*max(max(L0)),(1-S.c)*max(max(L0)));
  if S.bp
    L0 = L0-rconv2(L0,S.lo_filt);
  end
  en = upConv(im,z^2*S.lo_filt,'reflect1',[z z],[1 1],z*size(im))+L0;
% $$$   en = upConv(im,4*S.lo_filt,'reflect1',[2 2],[1 1],2*size(im))+L0;
