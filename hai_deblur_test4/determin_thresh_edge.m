function [threshold_r, threshold_s] =determin_thresh_edge(B, opt, thresh_percent)
    %set threshold to 90% data can be used in the code

    r0 = gradient_confidence_full(B,opt.window_size);
    sorted_r0 = sort(r0(:));
    nonzero_sorted_r0 = sorted_r0(sorted_r0>0);
    index_r = floor(size(nonzero_sorted_r0,1)*(1-thresh_percent));
    threshold_r = nonzero_sorted_r0(index_r);

    I_x0 = conv2(B, [-1,1;0,0], 'valid'); %vertical edges
    I_y0 = conv2(B, [-1,0;1,0], 'valid'); % horizontal edges         
    I_mag0 = sqrt(I_x0.^2+I_y0.^2);
    sorted_I_mag0 = sort(I_mag0(:));
    nonzero_sorted_I_mag0 = sorted_I_mag0(sorted_I_mag0>0);
    index_s = floor(size(nonzero_sorted_I_mag0,1)*(1-thresh_percent));
    threshold_s = nonzero_sorted_I_mag0(index_s);

end