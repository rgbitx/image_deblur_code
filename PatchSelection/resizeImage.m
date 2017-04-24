I = imread('image1.jpg');
Iscale = imresize(I,0.5);
imshow(Iscale)
imwrite(Iscale,'image1_1.png');