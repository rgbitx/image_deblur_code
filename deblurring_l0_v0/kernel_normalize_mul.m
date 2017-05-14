function K = kernel_normalize_mul(K,threshold)
if ~exist('threshold','var')
    threshold = 0;
end
% set elements below threshold to 0
[r,c,n] = size(K);
Kvec = reshape(K,r*c,[]);
max_Kvec = max(Kvec,[],1);
K(K<reshape(repmat(max_Kvec*threshold,r*c,1),[r,c,n])) = 0;
% sum is one
Kvec = reshape(K,r*c,[]);
K = K./reshape(repmat(sum(Kvec,1),r*c,1),[r,c,n]);% normalize

end

