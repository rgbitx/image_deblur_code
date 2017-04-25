function [k] = solve_kernel_irls(A,b,k_sz1,k_sz2,scla, lambda_kernel_smooth)
%% Note: This function is inefficient and will be replaced by our fast algorithm 

exp_a=2;
thr_0=0.0001; 

if ~exist('scla','var')
  scla=0.005;
end
A0=(A+A')/2; 
clear A;
opts = optimset('Display','off');
%k=quadprog(A0,-b,[],[],[],[],zeros(k_sz1*k_sz2,1),[],[],opts);
k = qpas(A0, -b,[],[],[],[],zeros(k_sz1*k_sz2,1));
%% fast implementation

% for itr=1:2%20 modified 2011-07-10 22:00:08
%   w=max(abs(k),thr_0).^(exp_a-2);
%   %k=quadprog(A0+scla*diag(w),-b,[],[],[],[],zeros(k_sz1*k_sz2,1),[],[],opts);
%    k = qpas(A0+scla*diag(w),-b,[],[],[],[],zeros(k_sz1*k_sz2,1));
%  %%
% end
% % % % % % % % % %  kernel=reshape(k,k_sz1,k_sz2);
% % % % % % % % % %  kernel = kernel-min(kernel(:));
% % % % % % % % % %  kernel = kernel/max(kernel(:));
% % % % % % % % % %  %%
% % % % % % % % % %  %for rmap test
% % % % % % % % % %  %k = smoothing_kernel(kernel, 1e-3, 1.5);
% % % % % % % % % %  k = smoothing_kernel(kernel, lambda_kernel_smooth, 1.5);
% % % % % % % % % % %   %for toy
% % % % % % % % % % %   if (max(k_sz1,k_sz2)>65)
% % % % % % % % % % %       k = smoothing_kernel(kernel, 1e-6, 3);
% % % % % % % % % % %   else
% % % % % % % % % % %       k = smoothing_kernel(kernel, 1e-6, 1.5);
% % % % % % % % % % %   end
 k(k<0) =0;
 k = k/sum(k(:));
k=reshape(k,k_sz1,k_sz2);
