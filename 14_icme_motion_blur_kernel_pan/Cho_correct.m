function kernel  = Cho_correct(kernel)
% pixel values are normalized into [0,1]
 ks = max(size(kernel));
 % to safely align the center of a kernel
 k1 = size(kernel)-1;
 kernel = padarray(kernel, k1);
 kernel = center_kernel(kernel);
 kernel = kernel(k1+1:end-k1, k1+1:end-k1);
end
 %-------------------------------------------
 % move the center of kernel
 % to the center of the kernel image
 function kernel = center_kernel(kernel)
 [kh kw] = size(kernel);
 [X Y] = meshgrid(1:kw, 1:kh);
 x_kernel_center = sum2(kernel .* X);
 y_kernel_center = sum2(kernel .* Y);
 x_img_center = (kw+1) / 2;
 y_img_center = (kh+1) / 2;
 x_shift = x_img_center - x_kernel_center;
 y_shift = y_img_center - y_kernel_center;
 kernel = warpimage(kernel, ...
 [1 0 -x_shift; 0 1 -y_shift]);
 end
 %-------------------------------------------
 function val = sum2(arr)
 val = sum(arr(:));
 end

 %-------------------------------------------
 % M: an inverse transform
 function warped = warpimage(im, M)
 [x,y]=meshgrid(1:size(im,2),1:size(im,1));
 coords=[x(:)'; y(:)'];
 homo_coords=[coords; ones(1,prod(size(im)))];
 warped_coords=M*homo_coords;
 xp=warped_coords(1,:);
 yp=warped_coords(2,:);

 warped = interp2(x,y,im,xp,yp,'linear');
 warped = reshape(warped,size(im));

 warped(isnan(warped))=0;
 end