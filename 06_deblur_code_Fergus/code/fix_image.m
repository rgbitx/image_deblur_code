function out = fix_image(in,reference);

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology
  
  SPACING = 0.05;
  
  %%% Assume input is image or patch of double format with +ve and -ve
  %%% values. Assume reference is +ve values, double, uint8 or uint16.    
    
  %%% histeq works bset when all histograms are [0:1]
  %%  so first, make reference image [0:1]

  
  ref_im = double(reference) / double(max(reference(:)));
  
  %% now take histogram of referecnce
  x_bins = [0:SPACING:1];
  hist_ref = hist(ref_im(:),x_bins);
    
  %%% now make input image in range 0 to 1 also
  m = min(in(:));
  in_shift = in - m;
  in_norm = in_shift / max(in_shift(:));
  
  %%% now use histeq
  out = histeq(in_norm,hist_ref);
  
    
    
