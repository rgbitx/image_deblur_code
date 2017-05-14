function K = kernel_initialize_mul(ks,kernel_type)
%output:
% K is ks X ks kernel
if ~exist('kernel_type','var')
    kernel_type = 'ver';
end
if ~exist('ks','var')
    ks = [3,3];
end
K = zeros(ks);
hks = floor(ks/2);
switch kernel_type
    case 'ver' %vertical
        K(hks(1)+1:hks(1)+2,hks(2)+1) = 1; 
    case 'hor' %horizontal
        K(hks(1)+1,hks(2)+1:hks(2)+2) = 1;
    case 'uniform'
        K(:) = 1;
    case 'delta'
        K(hks(1)+1,hks(2)+1) = 1;             
end% switch
sumK = sum(reshape(K,size(K,1)*size(K,2),[]),1);
K = K./reshape(repmat(sumK,size(K,1)*size(K,2),1),[size(K,1),size(K,2)]);% normalize

end

