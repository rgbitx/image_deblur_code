function S = L0Restoration_mul(B,kernel,lambda)
%%
% Image restoration with L0 prior
% The objective function: 
% S^* = argmin \sum_i^n 1/n ||I*k_i - B_i||^2 + lambda |\nabla I|_0
% qual to
% argmin \sum_i^n 1/n 1/lambda ||I*k_i - B_i||^2 + |\nabla I|_0
%% Input:
% @Im: Blurred image
% @kernel: blur kernel
% @lambda: weight for the L0 prior
% @kappa: Update ratio in the ADM
%% Output:
% @S: Latent image
%
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang,
%        Deblurring Text Images via L0-Regularized Intensity and Gradient
%        Prior, CVPR, 2014. 
%   [2] Li Xu, Cewu Lu, Yi Xu, and Jiaya Jia. Image smoothing via l0 gradient minimization.
%        ACM Trans. Graph., 30(6):174, 2011.
%
%   Author: Jinshan Pan (sdluran@gmail.com)
%   Date  : 05/18/2014

%% pad image
H = size(B,1);    W = size(B,2);
B = wrap_boundary_liu(B,opt_fft_size([H+size(kernel,1),W+size(kernel,2)]-1));
%% pre compuate
fx = [1, -1];
fy = [1; -1];
[N,M,D] = size(B);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
FB = B;
for ii = 1:D
    FB(:,:,ii) = fft2(B(:,:,ii));
end
%% 
n = D;
S = mean(B,3);
Den_KER = zeros(sizeI2D);
Normin1 = zeros(sizeI2D);
for ii = 1:D
    this_KER = psf2otf(kernel(:,:,ii),sizeI2D);
    Den_KER = Den_KER+abs(this_KER).^2/n;    
    Normin1 = Normin1+conj(this_KER).*FB(:,:,ii)/n;
end
%% 
betamax = 1e5;
beta = min(1,2*lambda);
beta_step = 2;
while beta < betamax
    %h-v
    Denormin  = Den_KER + beta*Denormin2;
    h = [diff(S,1,2), S(:,1,:)-S(:,end,:)];
    v = [diff(S,1,1); S(1,:,:)-S(end,:,:)];
    t = (h.^2+v.^2)<lambda/beta;    
    h(t)=0; 
    v(t)=0;
    % get S
    Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
    Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
    FS = (Normin1 + beta*fft2(Normin2))./Denormin;
    S = real(ifft2(FS));   
    % update beta
    beta = beta*beta_step;
end
S = S(1:H,1:W,:);

end