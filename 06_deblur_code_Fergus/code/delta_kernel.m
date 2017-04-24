function out = delta_kernel(s)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

  
  %%% if even size, make odd
  if (mod(s,2)==0)
    s = s + 1;
  end
  
  out = zeros(s);
  
  c = floor(s/2)+1;
  
  out(c,c)=1;
  
  
