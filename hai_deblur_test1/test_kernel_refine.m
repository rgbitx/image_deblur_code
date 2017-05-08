clear
close all

kernel = load('kernel_test_data.txt');

%%%%get epsilon_s
% threshold for epsilon_s
[width,height] = size(kernel);
k = kernel;

for it=1:5


threshD = 7 * norm(k,inf)/(2*width*height*5);
% threshD = norm(kernel,inf)/(2*width*it);

sorted_k = sort(k(:));
diff_sorted_k = diff(sorted_k);

all_index = find(diff_sorted_k>threshD);

if isempty(all_index)
   continue; 
end

epsilon_s = sorted_k(all_index(1)+1);

mask_s = heaviside_function(k,epsilon_s);
mask_s_ = ~mask_s;

s = mask_s .* k;
s_ = mask_s_ .* k;

% set up options for the kernel estimation
opts.use_fft = 1;
opts.lambda = 1;
opts.pcg_tol=1e-4;
opts.pcg_its = 1;
opts.s_ = s_;

k_prev = k;

k = pcg_kernel_irls_conv(k_prev, x1, y2, opts);

epsilon_k = (k-k_prev)/k;

figure, imagesc(s);

end

