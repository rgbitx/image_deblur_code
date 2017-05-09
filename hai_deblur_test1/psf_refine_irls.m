function psf = psf_refine_irls(blurred_x, blurred_y, latent_x, latent_y, kernel0, s_, weight, psf_size)

  %2013/5/31  %%%%%%%%%%%%%%%%%%%%%%%%%
    latent_xf = fft2(latent_x);
    latent_yf = fft2(latent_y);
    blurred_xf = fft2(blurred_x);
    blurred_yf = fft2(blurred_y);

    b_f = conj(latent_xf).* blurred_xf + conj(latent_yf).* blurred_yf;
    b = real(otf2psf(b_f, psf_size));

    p.m = conj(latent_xf).* latent_xf + conj(latent_yf).* latent_yf;

    p.img_size = size(blurred_xf);
    p.psf_size = psf_size;
    p.lambda = weight;
    p.s_ = s_;
    exp_a = 1;

%     psf = ones(psf_size) / prod(psf_size);
    psf = kernel0;

    for iter = 1:2
        psf_prev = psf;
        p.weightsl1 = p.lambda * (max(norm(psf_prev,1),1e-4) .^ (exp_a-2));                
        psf = conjgrad(psf, b, 20, 1e-5, @compute_Ax, p);

        psf(psf < max(psf(:))*0.05) = 0;
        psf = psf / sum(psf(:));
        
        figure,imagesc(psf);        
    end
end

function y = compute_Ax(x, p)
    x_f = psf2otf(x, p.img_size);
    y = otf2psf(p.m .* x_f, p.psf_size);
    y = y + p.weightsl1 * p.s_;
end
