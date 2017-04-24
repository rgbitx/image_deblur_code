clear
close all

%% initialization
I = imread('image1_1.png');
if size(I,3)==3
    Igray = rgb2gray(I);
end
y = im2double(Igray);
% y = Igray;

% y = shock_filter(y);

% figure;imshow(I,'InitialMagnification','fit');

ksize = 25;
psize = ksize * 16;

% calculate the gradients
% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];

gradient{1} = conv2(y, dx, 'valid');
gradient{2} = conv2(y, dy, 'valid');

fp = 0;
fpmax = fp;
fpmin = 1000;
pxmin=0;
pxmax=0;
pymin=0;
pymax=0;

%% calculate confidence
% tic
% calculate all confidences based on defined step
% pixAxis = offset + scale*(confAxis-1)
stepConf = floor(ksize);
[confs,offset,scale]=calconfidence(gradient,ksize,stepConf);
% toc


%% get the patch based on confidence
[rows,cols] = size(gradient{1});

stepPatch = floor(psize/3);
row=0;
col=0;

% tic
% go through all patches in a image
for row=psize:stepPatch:rows
    for col=psize:stepPatch:cols
              
        % pixAxis = offset + scale*(confAxis-1)
        xmin=row-psize+1;
        if xmin == 1
            cfxmin = 1;
        else
            cfxmin = floor((xmin - offset)/scale) + 1;
        end
        xmax=row;
        cfxmax = floor((xmax - offset)/scale);
        
        ymin=col-psize+1;
        if ymin == 1
            cfymin = 1;
        else
            cfymin = floor((ymin - offset)/scale) + 1;
        end
        ymax=col;
        cfymax = floor((ymax - offset)/scale);
        
        subConfidence = confs(cfxmin:cfxmax,cfymin:cfymax);
        
        fp = mean(subConfidence(:));
        
        if fp > fpmax
            fpmax = fp;
            
            % save the patch axis at max fp
            pxmin = xmin;
            pxmax = xmax;
            pymin = ymin;
            pymax = ymax;
            
            cfpxmin = cfxmin;
            cfpxmax = cfxmax;
            cfpymin = cfymin;
            cfpymax = cfymax;
            
        end
        
    end
end

assert(pxmin~=0 && pxmax~=0 && pymin~=0 && pymax~=0,...
       'pxmin,pxmax,pymin and pymax should be nonzero');

yp=y(pxmin:pxmax,pymin:pymax);
figure;imshow(yp,'InitialMagnification','fit');
% toc

%% calculate strong edge Is
threshold_r = 0.25;
threshold_s = 5e-3;

confPatch = confs(cfpxmin:cfpxmax,cfpymin:cfpymax);
M = heaviside(confPatch - threshold_r);

% unit zeros and ones for the mask
unit_zeros = zeros(stepConf);
unit_ones = ones(stepConf);

mask = ones






