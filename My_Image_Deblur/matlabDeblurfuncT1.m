LEN = 25;
THETA = 90;
PSF = fspecial('motion', LEN, THETA);
blurred = imread('t1.jpg');
wnr1 = deconvwnr(blurred, PSF, 0);
imshow(wnr1);