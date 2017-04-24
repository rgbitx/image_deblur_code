function [im,invKhat]=reconsEdge3(dx,dy,invKhat)

% Author: Yair Weiss
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% given dx dy reconstruct image
% here we do the convolutions using FFT2


im=zeros(size(dx));
[sx,sy]=size(dx);
mxsize=max(sx,sy);

if ~exist('invKhat')
  invK=invDel2(2*mxsize);
  invKhat=fft2(invK);
end

imX=conv2(dx,fliplr([0 1 -1]),'same');
imY=conv2(dy,flipud([0;1;-1]),'same');

imS=imX+imY;

imShat=fft2(imS,2*mxsize,2*mxsize);
im=real(ifft2(imShat.*invKhat));
im=im(mxsize+1:mxsize+sx,mxsize+1:mxsize+sy);




