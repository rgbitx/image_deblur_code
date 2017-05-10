function psf = adjust_psf_center(psf)

[X Y] = meshgrid(1:size(psf,2), 1:size(psf,1));
xc1 = sum2(psf .* X);
yc1 = sum2(psf .* Y);
xc2 = (size(psf,2)+1) / 2;
yc2 = (size(psf,1)+1) / 2;
xshift = round(xc2 - xc1);
yshift = round(yc2 - yc1);
psf = warpimage(psf, [1 0 -xshift; 0 1 -yshift]);

function val = sum2(arr)
val = sum(arr(:));