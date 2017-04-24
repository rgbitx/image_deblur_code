function [mask,mask_1] = compute_mask1(I,kernel_size)

  % 2016/06/01 %%%%%%%%%%%%%%%%%%%%%%%%%%
% compute mask

mask_1 = edge(I,'canny');
t = floor(kernel_size/1.5);

% t = floor(kernel_size/4);
se=strel('disk',t);
mask_1=imdilate(mask_1,se);

mask = ones(size(I))-mask_1;
