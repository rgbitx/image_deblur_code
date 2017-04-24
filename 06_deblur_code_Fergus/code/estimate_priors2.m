function [pi_x,gamma_x]=estimate_priors2(fname,imn,num,num_scales,type)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% Warning, this code is a terrible mess.

DEBUG = 1;  
SCALE_STEP = sqrt(2);  
GAMMA_CORRECTION = 1; 
INTENSITY_SCALING = 1/256;
IMAGE_OR_MAT = 'mat';
STEP_SIZE = 1;
PRESCALE = 1;
num_images = length(imn);
close all
MAX_IM_SIZE = 700;

if (isempty(imn))
  if strcmp(IMAGE_OR_MAT,'mat')
    imn{1} = '/data2/markov/blur/images/linear/sharp/image_0001.mat';
%    imn{1} = '../images/linear/whiteboard/image_0040.mat';   
    x = [-1:0.01:1];
  else
    imn{1} = '../images/aaron/cuba_sharp.jpg';
    %im{1} = imread('../images/linearblurry/image_0010.tif');
    x = [-200:1:200];
 end
 
 num_images = 1;

end
x = [-200:STEP_SIZE:200];

for b=1:num_scales
    
  scale = SCALE_STEP^(-(b-1));
  b_all = [];
 
  for a=1:num_images
    
    if strcmp(IMAGE_OR_MAT,'mat')
      load(imn{a});
      %x = [-1:0.01:1];
    else %% image
      im = imread(imn{a});
      %x = [-200:1:200];
    end 
    
    if (size(im,3)>1)
      im = double(rgb2gray(im));
    end
    
    %if (PRESCALE)
    %  im = imresize(im,PRESCALE,'bilinear');
    %end
    
    [imy,imx] = size(im);
    
    scale_factor = MAX_IM_SIZE / max([imx,imy]);
    
    im = imresize(im,scale_factor,'bilinear');
    
    im = im * INTENSITY_SCALING;
    
    if (GAMMA_CORRECTION~=1)
      im = (im.^GAMMA_CORRECTION)./(256^(GAMMA_CORRECTION-1));
    end
  
    if strcmp(type,'steer')
        
        [pyr,indices] = buildSpyr(im);
        
        b_x = spyrBand(pyr,indices,1,1);
        b_y = spyrBand(pyr,indices,1,2);
            
    elseif strcmp(type,'haar')
        
        b_x = conv2(im,[1 -1],'valid');
        b_y = conv2(im,[1 -1]','valid');
        
        if (scale~=1)
          b_x = imresize(b_x,scale,'bilinear');
          b_y = imresize(b_y,scale,'bilinear');          
        end
        
    else
        error foo
    end
 
    b_all = [ b_all , b_x(:)' , b_y(:)' ];
    
  end
  
  if DEBUG
    hista = hist(b_all,x);
    figure(b*2-1); plot(x,hista,'r'); hold on;
  end

  [mu,sigma,weight,log_likelihood]=GaussianMixtures1D(b_all,num);
    
  priors(b).pi = weight(1,:);
  priors(b).gamma = 1./ reshape(sigma(1,1,:),1,num);
  
  if DEBUG
    pdf_total_x = zeros(size(x));
    
    for c=1:num
      pdf_x = priors(b).pi(c) * normpdf(x,0,sqrt(1/priors(b).gamma(c)));
      pdf_total_x = pdf_total_x + pdf_x;
      figure(1); plot(x,pdf_x*STEP_SIZE*sum(hista),'b');   
    end     
    figure(b*2-1); plot(x,pdf_total_x*STEP_SIZE*sum(hista),'g');
  keyboard    
    figure(b*2);
    plot(x,log2((pdf_total_x)),'g');
    hold on;
    plot(x,log2((hista/(STEP_SIZE*sum(hista)))),'r');
  end
  
end

save(fname,'priors');
