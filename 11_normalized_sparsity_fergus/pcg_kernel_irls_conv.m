function k_out = pcg_kernel_irls_conv(k_init, X, Y, opts)
  
% 
% Use Iterative Re-weighted Least Squares to solve l_1 regularized kernel
% update with sum to 1 and nonnegativity constraints. The problem that is
% being minimized is:
%
% min 1/2\|Xk - Y\|^2 + \lambda \|k\|_1
% 
% Inputs:
% k_init = initial kernel, or scalar specifying size
% X = sharp image
% Y = blurry image  
% opts = options (see below)  
%
% Outputs:  
% k_out = output kernel
% 
% This version of the function uses spatial convolutions. Everything is maintained as 2D arrays
%
  
% Defaults
if nargin == 3
  opts.lambda = 0;  
  % PCG parameters
  opts.pcg_tol = 1e-8;
  opts.pcg_its = 100;
  fprintf('Input options not defined - really no reg/constraints on the kernel?\n');
end

% lambda = opts.lambda;
lambda = 1;
pcg_tol = opts.pcg_tol;
pcg_its = opts.pcg_its;
  
if (length(k_init(:)) == 1)
  k_init = zeros(k_init,k_init);
end;
  
% assume square kernel
ks = size(k_init,1);
ks2 = floor(ks/2); 
  
% precompute RHS
for i = 1:length(X)
  flipX{i} = fliplr(flipud(X{i}));
  % precompute X^T Y term on RHS (e = ks^2 length vector of all 1's)
  rhs{i} = conv2(flipX{i}, Y{i}, 'valid');
end;
  
tmp = zeros(size(rhs{1}));
for i = 1 : length(X)
  tmp = tmp + rhs{i};
end;
rhs = tmp;
  
k_out = k_init;
  
% Set exponent for regularization
exp_a = 1;

% outer loop
for iter = 1 : 1
  k_prev = k_out;
  % compute diagonal weights for IRLS
  weights_l1 = lambda .* (max(abs(k_prev), 0.0001) .^ (exp_a - 2));
  k_out = local_cg(k_prev, X, flipX, ks,weights_l1, rhs, pcg_tol, pcg_its);
end;
  
% local implementation of CG to solve the reweighted least squares problem
function k = local_cg(k, X, flipX, ks, weights_l1, rhs, tol, max_its)

Ak = pcg_kernel_core_irls_conv(k, X, flipX, ks,weights_l1);

r = rhs - Ak;

for iter = 1:max_its
  rho = (r(:)' * r(:));

  if (iter > 1)                      
    beta = rho / rho_1;
    p = r + beta*p;
  else
    p = r;
  end

  Ap = pcg_kernel_core_irls_conv(p, X, flipX, ks, weights_l1);

  q = Ap;
  alpha = rho / (p(:)' * q(:) );
  k = k + alpha * p;                    
  r = r - alpha*q;                    
  rho_1 = rho;
  
  if (rho < tol)
    break;
  end;
end;

  
    
    
  
  
   
