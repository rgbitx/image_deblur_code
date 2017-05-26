function S = L0regulation_xu(Im,kernel,lambda,epsilon)

if ~exist('epsilon','var')
    epsilon = 1;
end

%% pad image
H = size(Im,1);    W = size(Im,2);
Im = wrap_boundary_liu(Im, opt_fft_size([H W]+size(kernel)-1));

%%
S = Im;
fx = [1, -1];
fy = [1; -1];
[rows,cols] = size(Im);
sizeI2D = [rows,cols];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);

% Ix = conv2(Im, fx, 'valid'); %vertical edges
% Iy = conv2(Im, fy, 'valid');  
%%
F_KERNEL = psf2otf(kernel,sizeI2D);
F_D = abs(otfFx).^2 + abs(otfFy ).^2;
F_y = fft2(S);
Normin1 = conj(F_KERNEL).*F_y;
%%
for i=1:4
    beta = lambda/epsilon^2;
   for j = 1:1/epsilon
       Denormin = conj(F_KERNEL).*F_KERNEL + beta*F_D;
       
       % h-v subproblem
       h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
       v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
       
%        t=abs(h)<=epsilon;
%        h(t) = 0;
%        t=abs(v)<=epsilon;
%        h(t) = 0;
       t = (h.^2+v.^2)<epsilon^2; 
       h(t)=0; v(t)=0;
       clear t;
       
       % S subproblem       
       Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
       Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
%        Normin2 = conj(otfFx).*fft2(h) + conj(otfFy).*fft2(v);        
       numinator = Normin1+beta*fft2(Normin2);
%        numinator = Normin1+beta*Normin2;
       FS = numinator./Denormin;
       S = real(ifft2(FS));
   end
    epsilon = epsilon/2;
end

S = S(1:H, 1:W, :);



