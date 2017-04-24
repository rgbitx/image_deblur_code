function mpad = centerpad(m,padsize)
%input: m a matrix, padsize a 2x1 array.
%output: mpad, a zeropadded matrix with m in its "center"
%(that is, the convolution center is mapped to the convolution center
%of mpad). 
%size(mpad) = padsize;

[mx my] = size(m);
mpadx = padsize(1);
mpady = padsize(2);

% place m in the upper left minor of size [mx my] 
mpad = zeros(padsize);
mpad(1:mx,1:my) = m;

% now shift mpad so that m and mpad's convolution center will coincide. 

% find m and mpad's convolution centers.
mcenterx = floor(mx/2) + 1;
mcentery = floor(my/2) + 1;

mpadcenterx = floor(mpadx/2) + 1;
mpadcentery = floor(mpady/2) + 1;

% calculate the shift required 
% to map m's center to mpad's center.
mshiftx = mpadcenterx - mcenterx;
mshifty = mpadcentery - mcentery;

%perform the shift
mpad = circshift(mpad, [mshiftx mshifty]);
