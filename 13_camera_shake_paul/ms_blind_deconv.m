function [kernel,yorig,deblur] = ms_blind_deconv(fn, opts)

%
% Do multi-scale blind deconvolution given input file name and options
% structure opts. Returns a double deblurred image along with estimated
% kernel. Following the kernel estimation, a non-blind deconvolution is run.
%
% Copyright (2011): Dilip Krishnan, Rob Fergus, New York University.
%

if (isempty(fn))
   if (isempty(opts.blur))
      fprintf('No image provided in fn or opts.blur!!!\n');
      return;
   else
      y = opts.blur;
   end;
else
   y = im2double(imread(fn));
end;


% prescale the image if it's too big; kernel size is defined for the SCALED image
for k = 1:size(y, 3)
   y1(:, :, k) = imresize(y(:, :, k), opts.prescale, 'bilinear');
end;
y = y1;

% save off for non-blind deconvolution
yorig = y;

% gamma correct
y = y.^opts.gamma_correct;

% use a window to estimate kernel
if (~isempty(opts.kernel_est_win))
   w = opts.kernel_est_win;
else
   w(1) = 1; w(3) = size(y,1);
   w(2) = 1; w(4) = size(y,2);
end

if (size(y, 3) == 3)
   y = rgb2gray(y([w(1):w(3)], [w(2):w(4)], :));
end

b = zeros(opts.kernel_size);
bhs = floor(size(b, 1)/2);

% set kernel size for coarsest level - must be odd
minsize = max(3, 2*floor(((opts.kernel_size - 1)/16)) + 1);
fprintf('Kernel size at coarsest level is %d\n', minsize);

% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];

% l2 norm of gradient images
l2norm = 6;

resize_step = sqrt(2);
% determine number of scales
num_scales = 1;
tmp = minsize;
while(tmp < opts.kernel_size)
   ksize(num_scales) = tmp;
   num_scales = num_scales + 1;
   tmp = ceil(tmp * resize_step);
   if (mod(tmp, 2) == 0)
      tmp = tmp + 1;
   end;
end;
ksize(num_scales) = opts.kernel_size;


% blind deconvolution - multiscale processing
for s = 1:num_scales
   
   k1 = ksize(s); k2 = k1;
   if s == 1,
      k = init_kernel(ksize(1));
   else
      k = imresize(k, [k1 k2], 'bilinear');
   end
   k(k < 0) = 0;
   k = k./sum(k(:));
   
   % image size at this level
   r = floor(size(y, 1) * k1 / size(b, 1));
   c = floor(size(y, 2) * k2 / size(b, 2));
   if (s == num_scales)
      r = size(y, 1);
      c = size(y, 2);
   end
   
   % resize y according to the ratio of filter sizes
   ys = imresize(y, [r c], 'bilinear');
   g{1} = conv2(ys, dx, 'valid');
   g{2} = conv2(ys, dy, 'valid');
   %for gi = 1:2, rat(gi) = l2norm/norm(g{gi}(:)); end   
   
   if s == 1,
      ls = g;
   else
      ls{1} = imresize(ls{1}, size(g{1}), 'bilinear');
      ls{2} = imresize(ls{2}, size(g{2}), 'bilinear');
   end
   
   % upscale the estimated derivative image from previous level
   for gi = 1:2, 
   %   g{gi}    = rat(gi)*g{gi};       
      g{gi} = g{gi}/norm(g{gi}(:))*l2norm;
      ls{gi} = ls{gi}/norm(ls{gi}(:))*l2norm;
   end
   
   fprintf('Processing scale %d/%d; kernel size %dx%d; image size %dx%d\n', ...
      s, num_scales, k1, k2, r, c);
   
   if s == num_scales 
       opts.xk_iter = opts.final_scale_xk_iter; 
   end
   
   tic
   [ls, k] = ss_blind_deconv(g, ls, k, opts);
   toc
   % center the kernel
   [ls{1}, ls{2}, k] = center_kernel_separate(ls{1}, ls{2}, k);   
end
kernel = k;

padsize = bhs;
if opts.decon_kernel,
   
   if (opts.use_ycbcr)
      if (size(yorig, 3) == 3)
         ycbcr = rgb2ycbcr(yorig);
      else
         ycbcr = yorig;
      end
      opts.nb_alpha = 1;
   end
   
   if (opts.use_ycbcr == 1)
      ypad = padarray(ycbcr(:,:,1), [padsize padsize], 'replicate', 'both');
      for a = 1:4
         ypad = edgetaper(ypad, kernel);
      end;
      
      tmp = fast_deconv_bregman(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
      deblur(:, :, 1) = tmp(bhs + 1 : end - bhs, bhs + 1 : end - bhs);
      if (size(ycbcr, 3) == 3)
         deblur(:, :, 2:3) = ycbcr(:, :, 2:3);
         deblur = ycbcr2rgb(deblur);
      end;
   else
      for j = 1:3
         ypad = padarray(yorig(:, :, j), [1 1] * bhs, 'replicate', 'both');
%          for a = 1:4
%             ypad = edgetaper(ypad, kernel);
%          end;
         
         tmp = fast_deconv_bregman(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
         deblur(:, :, j) = tmp(bhs + 1 : end - bhs, bhs + 1 : end - bhs);
      end;
   end;
else
   deblur = [];
end


end

function [k] = init_kernel(minsize)
k = zeros(minsize, minsize);
k((minsize - 1)/2, (minsize - 1)/2:(minsize - 1)/2+1) = 1/2;
end