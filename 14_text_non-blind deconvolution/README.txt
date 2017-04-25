Non-blind image deconvolution



Run "non_blind_deconvolution_test.m" to try the examples included in this package.

===========================================
Parameters setting:

lambda_tv:     the weight for the Laplacian prior based deconvolution [1e-3,1e-2];
lambda_l0:     the weight for the L0 prior based deconvolution typically set as 2e-3, the best range is [1e-4, 2e-3].
weight_ring:   the larger values help suppress the ringing artifacts. weight_ring=0 imposes no suppression.

===========================================

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.).

This code and the algorithm are for non-comercial use only.


The implementation is based on the following paper:

Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang, "Deblurring Text Images via L0-Regularized Intensity and Gradient Prior", CVPR 2014


Contact: Jinshan Pan (sdluran@gmail.com)





