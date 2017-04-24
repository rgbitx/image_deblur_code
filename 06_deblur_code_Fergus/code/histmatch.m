function out = histmatch(in,reference)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology
  
  
  %%% in - double image btw. 0 and 1
  %%% ref - unit8 image btw. 0 and 255
  %%% out - uint8 image btw. 0 and 255  
  
  if ~strcmp(class(in),'double')
    error('Input image must be double');
  end

  if ~strcmp(class(reference),'uint8')
    error('Reference image must be UINT8');
  end  
   
  %if (size(in,3)~=3)
  %  error('Input image must be color');
  %end
  
  %if (size(reference,3)~=3)
  %  error('Reference image must be color');
  %end
    
  %%% Convert input images to grayscale
  if (size(in,3)~=1)
    gray_in = rgb2gray(in);
  else
    gray_in = in;
  end
  
  if (size(reference,3)~=1)
    gray_reference = rgb2gray(reference);
  else
    gray_reference = reference;
  end
  
  %%% Compute refernce histogram
  hist_reference = hist(double(gray_reference(:)),[0:255]);
  
  %%% equalize histograms
  [j,t] = histeq(gray_in,hist_reference);
  
  %% Now compute output image
  for a=1:size(in,3)
    %%% if change is small and round off won't have an effect
    %%% if change is large, all useful values will be compressed
    %%% to just a few levels, so this should be avoided
    %  out(:,:,a) = uint8(255*t(double(in(:,:,a))+1));

    %%% this version does everything as a double
    q = in(:,:,a);
    qm=interp1([0:255]/256,t,q(:));
    out(:,:,a) = uint8(256 * reshape(qm,size(in(:,:,a))));
    
  end
