function [structure] = structure_adaptive_map(im, theta, nIters)
% Decompose the input IMAGE into structure and texture parts using the
% Rudin-Osher-Fatemi method. The final output is a linear combination 
% of the decomposed texture and the structure parts.
%
% According to Wedel etal "An Improved Algorithm for TV-L1 optical flow"
%       equations (8)-(10)
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2009 $
if nargin == 1
    theta   = 1/8;
    nIters  = 100;
end;

% Rescale the input image to [-1 1]
IM   = scale_image(im, -1,1);

% Backup orginal images
im   = IM;

% stepsize
delta = 1.0./(4.0*theta+0.1);

for iIm = 1:size(im,3)
    
    % Initialize dual variable p to be 0
    p = zeros([size(im,1) size(im,2) 2]);
    
    % Gradient descend
    I = squeeze(IM(:,:,iIm));
    
    for iter = 1:nIters
        
        % Compute divergence
        div_p = imfilter(p(:,:,1), [-1 1 0], 'corr', 0)+ ...
            imfilter(p(:,:,2), [-1 1 0]', 'corr', 0);
        
        I_x = imfilter(I+theta.*div_p, [-1 1], 'replicate');
        
        I_y = imfilter(I+theta.*div_p, [-1 1]', 'replicate');
        
        % Update dual variable
        p(:,:,1) = p(:,:,1) + delta.*(I_x);
        p(:,:,2) = p(:,:,2) + delta.*(I_y);
        
        % Reproject to |p| <= 1 
        reprojection = max(1.0, sqrt(p(:,:,1).^2 + p(:,:,2).^2));
        p(:,:,1) = p(:,:,1)./reprojection;
        p(:,:,2) = p(:,:,2)./reprojection;
        
    end
    
    % compute divergence
    div_p = imfilter(p(:,:,1), [-1 1 0], 'corr', 0)+ ...
        imfilter(p(:,:,2), [-1 1 0]', 'corr', 0);
    
    % compute structure component
    IM(:,:,iIm) = I + theta.*div_p;
    
end;

if nargout == 1
    structure = squeeze(IM);
end;
end
%%
function imo = scale_image(im, vlow, vhigh, ilow, ihigh)
if nargin == 3
    ilow    = min(im(:));
    ihigh   = max(im(:));
end;

imo = (im-ilow)/(ihigh-ilow) * (vhigh-vlow) + vlow;
end
