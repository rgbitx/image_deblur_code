function lambda_scale = initial_weights(B)
%input:
% B: blurry images, in [0,255]
%output:
% Ms: mask for smooth region

if max(B(:)) < 2
    B = B*255;
end

threshold = 1;

n = size(B,3);
fsize = 5;
gfilter = fspecial('gaussian',[fsize,fsize],sqrt(2));
Bx = zeros(size(B));
By = zeros(size(B));
for ii = 1:n   
    % denoise using Gaussian filter
    B(:,:,ii) = imfilter(B(:,:,ii),gfilter,'same', 'replicate');
    Bx(:,:,ii) = imfilter(B(:,:,ii), [0 -1 1], 'same', 'replicate');
    By(:,:,ii) = imfilter(B(:,:,ii), [0;-1;1], 'same', 'replicate');    
end
% compute gradient
G = sqrt(Bx.^2+By.^2);
G = min(G,[],3);

lambda_scale = sum(G(:)<=threshold)/numel(G); %min_lambda

end