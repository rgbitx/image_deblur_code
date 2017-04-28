function [I] = deblurring_adm_aniso_1(B, k, lambda, alpha, lambda_1, kernel_size)

  % 2016/06/01 %%%%%%%%%%%%%%%%%%%%%%%%%%
% reference: Pan's code

beta = 1/lambda;

beta_min = 0.001;

[m n] = size(B);
% initialize with input or passed in initialization
I = B;

% make sure k is a odd-sized
if ((mod(size(k, 1), 2) ~= 1) | (mod(size(k, 2), 2) ~= 1))
    fprintf('Error - blur kernel k must be odd-sized.\n');
    return;
end;

[Nomin1, Denom1, Denom2] = computeDenominator(B, k);
Ix = [diff(I, 1, 2), I(:,1) - I(:,n)];
Iy = [diff(I, 1, 1); I(1,:) - I(m,:)];

%% Main loop
t = 0.01;
while beta > beta_min
    gamma = 1/(2*beta);
    Denom = Denom1 + gamma*Denom2;
    Wx = newton_w(Ix, 1/(beta*lambda), alpha, lambda_1, kernel_size);
    Wy = newton_w(Iy, 1/(beta*lambda), alpha, lambda_1, kernel_size);
    Wxx = [Wx(:,n) - Wx(:, 1), -diff(Wx,1,2)];
    Wxx = Wxx + [Wy(m,:) - Wy(1, :); -diff(Wy,1,1)];

    Fyout = (Nomin1 + gamma*fft2(Wxx))./Denom;
    I = real(ifft2(Fyout));

    Ix = [diff(I, 1, 2), I(:,1) - I(:,n)];
    Iy = [diff(I, 1, 1); I(1,:) - I(m,:)];
    beta = beta/2;
    t = t/1.5;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nomin1, Denom1, Denom2] = computeDenominator(y, k)
sizey = size(y);
otfk  = psf2otf(k, sizey);
Nomin1 = conj(otfk).*fft2(y);
Denom1 = abs(otfk).^2;
Denom2 = abs(psf2otf([1,-1],sizey)).^2 + abs(psf2otf([1;-1],sizey)).^2;
