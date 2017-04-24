clear

a = [1 1 1; 2 2 2; 3 3 3];
b =[1 0; 0 1; 1 1];

botf = psf2otf(b,[3,3]);
bfft = fft2(b,3,3);
afft = fft2(a);

yfft_oft = afft.*botf;
yfft_fft1 = bfft.*afft;
yfft_fft2 = afft.*bfft;
yifft_oft = real(ifft2(yfft_oft));
yifft_fft1 = real(ifft2(yfft_fft1));
yifft_fft2 = real(ifft2(yfft_fft2));

yconv1 = conv2(a,b,'valid');
yconv2 = conv2(a,b,'same');
yconv3 = conv2(b,a,'same');
