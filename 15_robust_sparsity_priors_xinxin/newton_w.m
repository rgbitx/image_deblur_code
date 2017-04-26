
function [w] = newton_w(v, beta, alpha, lambda_1, kernel_size)

  % 2016/06/01 %%%%%%%%%%%%%%%%%%%%%%%%%%

[mask,mask_1] = compute_mask1(v,kernel_size);
iterations = 4;

x = v;

for a=1:iterations
  fd = (alpha)*sign(x).*abs(x).^(alpha-1).*mask_1+beta*(x-v)+lambda_1*beta*mask.*x;
  fdd = alpha*(alpha-1)*abs(x).^(alpha-2).*mask_1+beta+lambda_1*beta*mask;
  x = x - fd./fdd;
end;

q = find(isnan(x));
x(q) = 0;
% check whether the zero solution is the better one
z = beta/2*v.^2;
f   = abs(x).^alpha + beta/2*(x-v).^2;
w = (f<z).*x;
