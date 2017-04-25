function r = edge_map(input,win_size)
% The function is to obtain selective edge map based on 
% Two-Phase Kernel Estimation for Robust Motion Deblurring by L. Xu and J. Jia ECCV 2010
% It is trying to rule out small scale edge, and returns zeros for small
% scale edge or smooth region, ones otherwise.
% 
% Param:
% input: the blurred image
% win_size: the size of the small edge been ruled out, odd number
% thresh: the threshold for r map

% if ~isdouble(input)
%     input = im2double(input);
% end
%input = double(uint8(input));

[Bx,By] = gradient(input);
Bx_abs = abs(Bx);
By_abs = abs(By);
mask = ones(win_size,win_size);

Bx_mask = conv2(Bx,mask,'same');
By_mask = conv2(By,mask,'same');
Bx_abs_mask = conv2(Bx_abs,mask,'same');
By_abs_mask = conv2(By_abs,mask,'same');
 r = abs(Bx_mask+By_mask)./(Bx_abs_mask+By_abs_mask+0.5);
%r = (abs(Bx_mask)+abs(By_mask))./(Bx_abs_mask+By_abs_mask+0.5);

% [m,n] = size(input);
% output = zeros(m,n);
% output(r>thresh) = 1;
% ratio = sum(output(:))/m/n;
