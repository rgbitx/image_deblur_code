function [yorig, deblur, kernel, opts] = ms_blind_deconv(fn, opts)

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
  if (size(y, 3) == 3)
    y = rgb2gray(y([w(1):w(3)], [w(2):w(4)], :));
  end;
else
  if (size(y, 3) == 3)
    y = rgb2gray(y);
  end;
end;

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
  if (s == 1)
    % at coarsest level, initialize kernel
    ks{s} = init_kernel(ksize(1));
    k1 = ksize(1);
    k2 = k1; % always square kernel assumed
  else
    % upsample kernel from previous level to next finer level
    k1 = ksize(s);
    k2 = k1; % always square kernel assumed
    
    % resize kernel from previous level
    tmp = ks{s-1};
    tmp(tmp<0) = 0;
    tmp = tmp/sum(tmp(:));
    ks{s} = imresize(tmp, [k1 k2], 'bilinear');
    % bilinear interpolation not guaranteed to sum to 1 - so renormalizes
    ks{s}(ks{s} < 0) = 0;
    sumk = sum(ks{s}(:));
    ks{s} = ks{s}./sumk;
  end;
  
  % image size at this level
  r = floor(size(y, 1) * k1 / size(b, 1));
  c = floor(size(y, 2) * k2 / size(b, 2));
  
  if (s == num_scales)
    r = size(y, 1);
    c = size(y, 2);
  end;
  
  fprintf('Processing scale %d/%d; kernel size %dx%d; image size %dx%d\n', ...
            s, num_scales, k1, k2, r, c);
  
    
  % resize y according to the ratio of filter sizes
  ys = imresize(y, [r c], 'bilinear');
  yx = conv2(ys, dx, 'valid'); 
  yy = conv2(ys, dy, 'valid'); 
  
  c = min(size(yx, 2), size(yy, 2));
  r = min(size(yx, 1), size(yy, 1));
  
  g = [yx yy];
  
  % normalize to have l2 norm of a certain size
  tmp1 = g(:, 1:c);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  g(:, 1:c) = tmp1;
  tmp1 = g(:, c+1:end);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  g(:, c+1:end) = tmp1;
  
  if (s == 1)
    ls{s} = g;
  else
    if (error_flag ~= 0)
      ls{s} = g;
    else
      % upscale the estimated derivative image from previous level
      c1 = (size(ls{s - 1}, 2)) / 2;
      tmp1 = ls{s - 1}(:, 1:c1);
      tmp1_up = imresize(tmp1, [r c], 'bilinear');
      tmp2 = ls{s - 1}(:, c1 + 1 : end);
      tmp2_up = imresize(tmp2, [r c], 'bilinear');
      ls{s} = [tmp1_up tmp2_up];
    end;
  end;

  tmp1 = ls{s}(:, 1:c);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  ls{s}(:, 1:c) = tmp1;
  tmp1 = ls{s}(:, c+1:end);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  ls{s}(:, c+1:end) = tmp1;
  
  % call kernel estimation for this scale
  opts.lambda{s} = opts.min_lambda;
  
  [ls{s} ks{s} error_flag] = ss_blind_deconv(g, ls{s}, ks{s}, opts.lambda{s}, ...
                                               opts.delta, opts.x_in_iter, opts.x_out_iter, ...
                                               opts.xk_iter, opts.k_reg_wt);

  if (error_flag < 0)
    ks{s}(:) = 0;
    ks{s}(ceil(size(ks{s}, 1)/2), ceil(size(ks{s}, 2)/2)) = 1;
    fprintf('Bad error - just set output to delta kernel and return\n');
  end;
  
   % center the kernel
  c1 = (size(ls{s}, 2)) / 2;
  tmp1 = ls{s}(:, 1:c1);
  tmp2 = ls{s}(:, c1 + 1 : end);
  [tmp1_shifted tmp2_shifted ks{s}] = center_kernel_separate(tmp1, tmp2, ks{s});
  ls{s} = [tmp1_shifted tmp2_shifted];
  
  % set elements below threshold to 0
  if (s == num_scales)
    kernel = ks{s};
    kernel(kernel(:) < opts.k_thresh * max(kernel(:))) = 0;
    kernel = kernel / sum(kernel(:));
  end;
end;

padsize = bhs;
dx = [1 -1];
dy = dx';
  
if (opts.use_ycbcr)
  if (size(yorig, 3) == 3)
    ycbcr = rgb2ycbcr(yorig);
  else
    ycbcr = yorig;
  end;
  opts.nb_alpha = 1; 
end;

if (opts.use_ycbcr == 1)
    ypad = padarray(ycbcr(:,:,1), [padsize padsize], 'replicate', 'both');
%     for a = 1:4 
%       ypad = edgetaper(ypad, kernel); 
%     end;   
    
    tmp = fast_deconv_bregman(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
    deblur(:, :, 1) = tmp(bhs + 1 : end - bhs, bhs + 1 : end - bhs);
    if (size(ycbcr, 3) == 3)
      deblur(:, :, 2:3) = ycbcr(:, :, 2:3);
      deblur = ycbcr2rgb(deblur);
    end;
else
  for j = 1:3
    ypad = padarray(yorig(:, :, j), [1 1] * bhs, 'replicate', 'both');
%     for a = 1:4 
%       ypad = edgetaper(ypad, kernel); 
%     end;
    
    tmp = fast_deconv_bregman(ypad, kernel, opts.nb_lambda, opts.nb_alpha);
    deblur(:, :, j) = tmp(bhs + 1 : end - bhs, bhs + 1 : end - bhs);
  end;
end;

figure; imagesc([uint8(255*yorig) uint8(255*deblur)]); title(['Blurred/' ...
                    'deblurred']);
figure; imagesc(kernel); colormap gray; title('Kernel');

function [k] = init_kernel(minsize)
  k = zeros(minsize, minsize);
  k((minsize - 1)/2, (minsize - 1)/2:(minsize - 1)/2+1) = 1/2;
