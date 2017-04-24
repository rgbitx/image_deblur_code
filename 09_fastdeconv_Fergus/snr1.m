function [x,x2] = snr1(sig, hk, ref)
% x = snr(sig, ref)
% snr -- Compute Signal-to-Noise Ratio for images
%
% Usage:
%       x = snr(sig, ref)  -- 1st time or
%       x = snr(sig)       -- afterwards
%
% Input:
%       sig         Modified image
%       ref         Reference image
%  
% Output:
%       x           SNR value
%       x2          SSD value

persistent ref_save;

if nargin == 3; ref_save = ref; end;

%% just use central portion
hk = hk*4;
ref_save2 = ref_save(hk+1:end-hk,hk+1:end-hk);
sig2 = sig(hk+1:end-hk,hk+1:end-hk);

mse = mean((ref_save2(:)-sig2(:)).^2);

dv = var(ref_save2(:),1);
x = 10*log10(dv/mse);

if nargout==2
  x2 = mse*prod(size(sig));
end
 
