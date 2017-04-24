function [pxmin,pxmax,pymin,pymax] = patchSelection(img, ksize)

I = imread(img);
Igray = rgb2gray(I);
y = im2double(Igray);

% figure;imshow(I,'InitialMagnification','fit');

% ksize = 25;
psize = ksize * 15;

% calculate the gradients
% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];

gradient{1} = conv2(y, dx, 'valid');
gradient{2} = conv2(y, dy, 'valid');

fp = 0;
fpmax = fp;
pxmin=0;
pxmax=0;
pymin=0;
pymax=0;

[rows,cols] = size(gradient{1});
% tic
% go through all patches in a image
for row=psize:psize:rows
    for col=psize:psize:cols
       
        xmin=row-psize+1;
        xmax=row;
        ymin=col-psize+1;
        ymax=col;
        
        subgradient = {gradient{1}(xmin:xmax,ymin:ymax),...
            gradient{2}(xmin:xmax,ymin:ymax)};
        
        fp = calconfidence(subgradient,ksize);
        
        if fp > fpmax
            fpmax = fp;
            
            pxmin = xmin;
            pxmax = xmax;
            pymin = ymin;
            pymax = ymax;
            
        end
        
    end
end

% yf=y(pxmin:pxmax,pymin:pymax);
% figure;imshow(yf,'InitialMagnification','fit');
% time=toc

