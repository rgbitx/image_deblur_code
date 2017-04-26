function r = gradient_confidence_full(B,d)
%% Compate rmap
% Reference:  
% Li Xu, Jiaya Jia, 
% Two-Phase Kernel Estimation for Robust Motion Deblurring,  ECCV 2010.

if (size(B,3) ==3)
    B0 = rgb2gray(B);
else
    B0 = B;
end
border = (d-1)/2;
B0 = double(B0);
y_padded = padarray(B0, [border border], 'symmetric', 'both');
B_X = conv2(y_padded,  [-1,1;0,0], 'same'); % vertical edges
B_Y = conv2(y_padded, [-1,0;1,0], 'same'); % horizontal edges
B_mag = sqrt(B_X.^2 + B_Y.^2); % gradient magnitude
B_X_block_sum = boxfilter(B_X, border);
B_Y_block_sum = boxfilter(B_Y, border);
B_mag_sum = boxfilter(B_mag, border);
r = sqrt(B_X_block_sum.^2 + B_Y_block_sum.^2)./(B_mag_sum + 0.5);
%r  = r (border+1:end-border, border+1:end-border);
r  = r (1:end-2*border, 1:end-2*border);
end
%% subfunction
function imDst = boxfilter(imSrc, r)

%   BOXFILTER   O(1) time box filtering using cumulative sum

[hei, wid] = size(imSrc);
imDst = zeros(size(imSrc));

%cumulative sum over Y axis
imCum = cumsum(imSrc, 1);
%difference over Y axis
imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);

%cumulative sum over X axis
imCum = cumsum(imDst, 2);
%difference over Y axis
imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
end


