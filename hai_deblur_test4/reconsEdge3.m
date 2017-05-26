function im = reconsEdge3(dx,dy)

% Author: Yair Weiss
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% Modified by Oliver Whyte for CVPR 2010 paper:
% "Non-uniform Deblurring for Shaken Images"
% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce

% given dx dy reconstruct image
% here we do the convolutions using FFT2

% if only one input, assume left half is gx, right half is gy
if nargin < 2, dy = dx(:,end/2+1:end,:); dx = dx(:,1:end/2,:); end


im=zeros(size(dx));
[sx,sy]=size(dx);
mxsize=max(sx,sy);

invK=invDel2(2*mxsize);
invKhat=fft2(invK);

imX=conv2(dx,fliplr([0 1 -1]),'same');
imY=conv2(dy,flipud([0;1;-1]),'same');

imS=imX+imY;

imShat=fft2(imS,2*mxsize,2*mxsize);
im=real(ifft2(imShat.*invKhat));
im=im(mxsize+1:mxsize+sx,mxsize+1:mxsize+sy);




function [invK]=invDel2(isize);

 K=zeros(isize);
 K(isize/2,isize/2)=-4;
 K(isize/2+1,isize/2)=1;
 K(isize/2,isize/2+1)=1;
 K(isize/2-1,isize/2)=1;
 K(isize/2,isize/2-1)=1;
 
 Khat=fft2(K);
 I=find(Khat==0);
 Khat(I)=1;
 invKhat=1./Khat;
 invKhat(I)=0;
 invK=ifft2(invKhat);
 invK=-real(invK);
 invK=conv2(invK,[1 0 0;0 0 0;0 0 0],'same');% shift by one
 
