function [g] = fast_deconv_bregman(f, k, lambda, alpha)

%
% fast solver for the non-blind deconvolution problem: min_g \lambda/2 |g \oplus k
% - f|^2. We use a splitting trick as
% follows: introduce a (vector) variable w, and rewrite the original
% problem as: min_{g,w,b} \lambda/2 |g \oplus k - g|^2 + \beta/2 |w -
% \nabla g - b|^2, and then we use alternations on g, w
% and b to update each one in turn. b is the Bregman variable. beta is
% fixed. An alternative is to use continuation but then we need to set a
% beta regime. Based on the NIPS 2009 paper of Krishnan and Fergus "Fast
% Image Deconvolution using Hyper-Laplacian Priors"
%
  
beta = 400;

initer_max = 1;
outiter_max = 5;

[m n] = size(f); 

% initialize
g = f; 

% make sure k is odd-sized
if ((mod(size(k, 1), 2) ~= 1) | (mod(size(k, 2), 2) ~= 1))
  fprintf('Error - blur kernel k must be odd-sized.\n');
  return;
end;
ks = floor((size(k, 1)-1)/2);

dx = [1 -1];
dy = dx';
dxt = fliplr(flipud(dx));
dyt = fliplr(flipud(dy));

[Ktf, KtK, DtD, Fdx, Fdy] = computeConstants(f, k, dx, dy);

gx = conv2(g, dx, 'valid');
gy = conv2(g, dy, 'valid');

fx = conv2(f, dx, 'valid');
fy = conv2(f, dy, 'valid');

ks = size(k, 1);
ks2 = floor(ks / 2);

% store some of the statistics
lcost = [];
pcost = [];
outiter = 0;

bx = zeros(size(gx));
by = zeros(size(gy));
wx = gx;
wy = gy;

totiter = 1;
gk = conv2(g, k, 'same');
      
lcost(totiter) = (lambda / 2) * norm(gk(:) - f(:))^2;
pcost(totiter) = sum((abs(gx(:)) .^ alpha));
pcost(totiter) = pcost(totiter) + sum((abs(gy(:)) .^ alpha));

for outiter = 1:outiter_max
  fprintf('Outer iteration %d\n', outiter);
  initer = 0;

  for initer = 1:initer_max
    totiter = totiter + 1;
      
    if (alpha == 1)
      tmpx = beta * (gx + bx);
      betax = beta;
      tmpx = tmpx ./ betax;
      
      tmpy = beta * (gy + by);
      betay = beta;
      tmpy = tmpy ./ betay;
      betay = betay;
      wx = max(abs(tmpx) - 1 ./ betax, 0) .* sign(tmpx);
      wy = max(abs(tmpy) - 1 ./ betay, 0) .* sign(tmpy);
    else
      wx = solve_image_bregman(gx + bx, beta, alpha); 
      wy = solve_image_bregman(gy + by, beta, alpha);
    end
    
    bx = bx - wx + gx;
    by = by - wy + gy;
      
    wx1 = conv2(wx - bx, dxt, 'full');
    wy1 = conv2(wy - by, dyt, 'full');
    tmp = zeros(size(g));
    
    gprev = g;
    gxprev = gx;
    gyprev = gy;
    
    num = lambda * Ktf + beta * fft2(wx1 + wy1);
    denom = lambda * KtK + beta * DtD;
    Fg = num ./ denom;
    g = real(ifft2(Fg));
    
    gx = conv2(g, dx, 'valid');
    gy = conv2(g, dy, 'valid');
    gk = conv2(g, k, 'same');
    lcost(totiter) = (lambda / 2) * norm(gk(:) - f(:))^2;
    pcost(totiter) = sum((abs(gx(:)) .^ alpha));
    pcost(totiter) = pcost(totiter) + sum((abs(gy(:)) .^ alpha));
  end
end

function [Ktf, KtK, DtD, Fdx, Fdy] = computeConstants(f, k, dx, dy)

sizef = size(f);
otfk  = psf2otf(k, sizef); 
Ktf = conj(otfk) .* fft2(f);
KtK = abs(otfk) .^ 2; 
Fdx = abs(psf2otf(dx, sizef)).^2;
Fdy = abs(psf2otf(dy, sizef)).^2;
DtD = Fdx + Fdy;
