function a = rgb2gray_rob(r,g,b)

% Author: Rob Fergus, adapted from Mathwords version
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

%%% Rob's version sets output to 255 if ANY of the channels are saturated
%%% Assumes a single input  
SATURATION_LEVEL = 250;  

%RGB2GRAY Convert RGB image or colormap to grayscale.
%   RGB2GRAY converts RGB images to grayscale by eliminating the
%   hue and saturation information while retaining the
%   luminance.
%
%   I = RGB2GRAY(RGB) converts the truecolor image RGB to the
%   grayscale intensity image I.
%
%   NEWMAP = RGB2GRAY(MAP) returns a grayscale colormap
%   equivalent to MAP.
%
%   Class Support
%   -------------
%   If the input is an RGB image, it can be of class uint8, 
%   uint16 or double; the output image I is of the same class 
%   as the input image. If the input is a colormap, the input 
%   and output colormaps are both of class double.
%
%   See also IND2GRAY, NTSC2RGB, RGB2IND, RGB2NTSC.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.20 $  $Date: 2002/03/15 15:29:04 $

%% get ind of pixels who exceed level in ANY channel
sat_mask = find((r(:,:,1)>SATURATION_LEVEL) | (r(:,:,2)>SATURATION_LEVEL) | (r(:,:,3)>SATURATION_LEVEL));

if nargin==0,
  error('Need input arguments.');
end
threeD = (ndims(r)==3); % Determine if input includes a 3-D array 

if nargin==1,
  if threeD,
    rgb = reshape(r(:),size(r,1)*size(r,2),3);
    a = zeros([size(r,1), size(r,2)]);
  else % Colormap
    error foo
  end

else
  error('Invalid input arguments.');
end

%%% find saturated pixels
%lowr = find(rgb(:,

% TODO

T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);

if isa(rgb, 'uint8')  
    a = uint8(reshape(double(rgb)*T(1,:)', size(a)));
elseif isa(rgb, 'uint16')
    a = uint16(reshape(double(rgb)*T(1,:)', size(a)));
elseif isa(rgb, 'double')    
    a = reshape(rgb*T(1,:)', size(a));
    a = min(max(a,0),1);
end

if ((nargin==1) & (~threeD)),    % rgb2gray(MAP)
    if ~isa(a, 'double')
        a = im2double(a);
    end
    a = [a,a,a];
end

%%% set sat pixels to 255
a(sat_mask) = 255;
