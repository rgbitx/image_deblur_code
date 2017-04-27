function [k] = solve_kernel_irls(A,b,k_sz1,k_sz2,scla, lambda_kernel_smooth)
%% Solve kernel by using IRLS method

exp_a=0.5;
thr_0=0.0001; 

if ~exist('scla','var')
  scla=0.005;
end
A0=(A+A')/2; 
clear A;
opts = optimset('Display','off');
k=quadprog(A0,-b,[],[],[],[],zeros(k_sz1*k_sz2,1),[],[],opts);
% k = qpas(A0, -b,[],[],[],[],zeros(k_sz1*k_sz2,1));

for itr=1:2 %20 modified 2011-07-10 22:00:08
  w=max(abs(k),thr_0).^(exp_a-2);
  k=quadprog(A0+scla*diag(w),-b,[],[],[],[],zeros(k_sz1*k_sz2,1),[],[],opts);
%    k = qpas(A0+scla*diag(w),-b,[],[],[],[],zeros(k_sz1*k_sz2,1));
 %%
end
 kernel=reshape(k,k_sz1,k_sz2);
 kernel = kernel-min(kernel(:));
 kernel = kernel/max(kernel(:));
 %%
 %for rmap test
 %k = smoothing_kernel(kernel, 1e-3, 1.5);
 if lambda_kernel_smooth~=0
     k = smoothing_kernel(kernel, lambda_kernel_smooth, 1.5);
 end

 k(k<0) =0;
 k = k/sum(k(:));
k=reshape(k,k_sz1,k_sz2);
