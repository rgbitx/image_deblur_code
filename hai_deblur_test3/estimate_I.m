function latent = estimate_I(Bp, kernel, weight)

    %----------------------------------------------------------------------
    % these values can be pre-computed at the beginning of each level
%     blurred_f = fft2(blurred);
%     dx_f = psf2otf([1 -1 0], size(blurred));
%     dy_f = psf2otf([1;-1;0], size(blurred));
%     blurred_xf = dx_f .* blurred_f; %% FFT (Bx)
%     blurred_yf = dy_f .* blurred_f; %% FFT (By)
    
    latent_size = size(Bp);

    kernel_f = psf2otf(kernel,latent_size);
    blurred_f = fft2(Bp);
    lambda = weight;
    
    dx = [-1 1; 0 0];
    dy = [-1 0; 1 0];
    
    dx_f = psf2otf(dx,latent_size);
    dy_f = psf2otf(dy,latent_size);
    
    % compute b = sum_i w_i latent_i * blurred_i
    b_f = conj(kernel_f)  .* blurred_f;
    b = real(otf2psf(b_f, latent_size));

    p.m = conj(kernel_f)  .* kernel_f ...
        + lambda*(conj(dx_f).*dx_f + conj(dy_f).*dy_f);
    %p.img_size = size(blurred);
    p.img_size = latent_size;
    p.psf_size = latent_size;
    p.lambda = weight;

    latent = ones(latent_size) / prod(latent_size);
    latent = conjgrad(latent, b, 20, 1e-5, @compute_Ax, p);
    
%     latent(latent < max(latent(:))*0.05) = 0;
%     latent = latent / sum(latent(:));    
end

function y = compute_Ax(x, p)
    x_f = psf2otf(x, p.img_size);
    y = otf2psf(p.m .* x_f, p.psf_size);
end
