
function [G]=shock_filter(I, maskSize, iterations, dt)

close all;

if (~exist('I','var'))
	I=rgb2gray(imread('C:\Data\leaves_small.jpg'));
end

I=double(I);

if (~exist('dt','var'))
	dt=0.25;
end

if (~exist('maskSize','var'))
	maskSize=9;
end

if (~exist('iterations','var'))
	iterations=10;
end

g=fspecial('gaussian',maskSize);
G=I;
%G=conv2(I,g,'same');
tic;
for i=1:iterations

    G=conv2(G,g,'same');
    [gx,gy]=gradient(G);
    normxy=sqrt(gx.*gx+gy.*gy);
    
    s = -sign(del2(G));    
    G=G+dt.*s.*normxy;
    
end
toc
figure,imshow(I);
figure,imshow(G);

end
