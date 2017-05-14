function [B_org,B] = preparedata(img_name)
B_org = imread(img_name);

[rr,cc] = size(B_org);
B = zeros(rr,cc,n,'uint8');

B(:,:,ii) = rgb2gray(B_org);


end