img = imread('image3.png');

y = im2double(img);

% convert to gray
if size(y,3)==3
    y = rgb2gray(y);
end

% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];

% calculate the gradient of image 
grad{1} = conv2(y,dx,'valid');
grad{2} = conv2(y,dy,'valid');



