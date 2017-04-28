This package contains an implementation of the image deblurring algorithm described in the paper: 
Jinshan Pan, Deqing Sun, Hanspeter Pfister, and Ming-Hsuan Yang, "Blind Image Deblurring Using Dark Channel Prior", CVPR 2016. 

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) 
in an academic publication.

For algorithmic details, please refer to our paper.

----------------
How to use
----------------
The code is tested in MATLAB 2013b(64bit) under the MS Windows 7 64bit version with an Intel Core i7-4790 CPU and 28 GB RAM.

1. unpack the package
2. include code/subdirectory in your Matlab path
3. Run "demo_deblurring.m" to try the examples included in this package.
----------------
User specified parameter:
----------------
There are a few parameters need to be specified by users.
---------------
Kernel estimation part:
---------------
'gamma_correct': gamma correction for the input image (typically set as 1 and 2.2)
'kernel_size':   the size of blur kernel
'lambda_dark':   the weight for the L0 regularization on dark channel prior (typically set as 4e-3)
'lambda_grad':   the weight for the L0 regularization on gradient (typically set as 4e-3)
---------------
Non-blind deconvolution part:
Our algorithm is mainly designed for blur kernel esitmation part. For the non-blind deblurring, we use existing methods (e.g., [1, 2, 3, 4]) to genetate the final results for fair comparisons.
[1] J. Pan, Z. Hu, Z. Su, and M.-H. Yang, Deblurring Text Images via L0-Regularized Intensity and Gradient Prior. CVPR 2014.
[2] S. Cho, J. Wang, and S. Lee, Handling outliers in non-blind image deconvolution. ICCV 2011.
[3] O. Whyte, J. Sivic, and A. Zisserman, Deblurring shaken and partially saturated images. ICCV Workshops 2011
[4] Z. Hu, S. Cho, J. Wang, and M.-H. Yang, Deblurring lowlight images with light streaks. CVPR 2014.
---------------

----------------
IMPORTANT NOTE 
----------------
Note that the algorithm sometimes may converge to an incorrect result. When you obtain such an incorrect result, please re-try to deblur with a slightly changed parameters (e.g., using large blur kernel sizes or gamma correction (2.2)). Should you have any questions regarding this code and the corresponding results, please contact Jinshan Pan (sdluran@gmail.com), Deqing Sun (deqings@nvidia.com).
