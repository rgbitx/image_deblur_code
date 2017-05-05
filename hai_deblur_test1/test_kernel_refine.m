clear

kernel = load('kernel_test_data.txt');

%%%%get epsilon_s
% threshold for epsilon_s
[width,height] = size(kernel);
threshD = 7 * norm(kernel,inf)/(2*width*height*5);

sorted_k = sort(kernel(:));
diff_sorted_k = diff(sorted_k);

all_index = find(diff_sorted_k>threshD);

assert(~isempty(all_index),'index should be non empty!');

epsilon_s = sorted_k(all_index(1)+1);

% mask_s = epsilon_s*ones(width,height);

S = kernel(kernel>epsilon_s);