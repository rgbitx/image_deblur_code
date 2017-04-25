function out = pcg_kernel_core_irls_conv(k, X, flipX, ks, weights_l1)

% This function applies the left hand side of the IRLS system to the
% kernel x. Uses conv2's. 
  
% first term: X'*X*k (quadratic term)
out_l2 = zeros(size(k));

for i = 1:length(X)
  tmp1 = conv2(X{i}, k, 'valid');
  tmp2 = conv2(flipX{i}, tmp1, 'valid');
  out_l2 = out_l2 + tmp2;
end;
  
% second term: L1 regularization
out_l1 = weights_l1 .* k;
  
out = out_l2 + out_l1;
