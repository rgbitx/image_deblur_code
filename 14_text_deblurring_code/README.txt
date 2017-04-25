Title:      Matlab code for "Deblurring Text Images via L0-Regularized Intensity and Gradient Prior"  
Author:     Jinshan Pan (sdluran@gmail.com), Zhe Hu (zhu@ucmerced.edu), Zhixun Su (zxsu@dlut.edu.cn), Ming-Hsuan Yang (mhyang@ucmerced.edu).  
Version:    1.0 
Copyright:  2014, Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang


Matlab code for "Deblurring Text Images via L0-Regularized Intensity and Gradient Prior"
==========================================================

This package contains an implementation of the image deblurring algorithm described in the paper: 
Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang, "Deblurring Text Images via L0-Regularized Intensity and Gradient Prior", CVPR 2014.

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) 
in an academic publication.

For algorithmic details, please refer to our paper.

----------------
How to use
----------------
The code is tested in MATLAB 2011b(32bit) under the MS Windows 7 64bit version with an Intel Xeon CPU@2.53GHz and 12GB RAM.

1. unpack the package
2. include code/subdirectory in your Matlab path
3. Run "demo_text_deblurring.m" to try the examples included in this package.
----------------
User specified parameter:
----------------
There are a few parameters need to be specified by users.
---------------
Kernel estimation part:
---------------
'gamma_correct': gamma correction for the input image (typically set as 1 and 2.2)
'kernel_size':   the size of blur kernel
'lambda_pixel':  the weight for the L0 regularization on intensity (typically set as 4e-3)
'lambda_grad':   the weight for the L0 regularization on gradient (typically set as 4e-3)
---------------
Non-blind deconvolution part:
---------------
'lambda_tv':     the weight for the Laplacian prior based deconvolution [1e-3,1e-2];
'lambda_l0':     the weight for the L0 prior based deconvolution typically set as 2e-3, the best range is [1e-4, 2e-3].
'weight_ring':   the larger values help suppress the ringing artifacts. weight_ring=0 imposes no suppression.

----------------
IMPORTANT NOTE 
----------------
Note that the algorithm sometimes may converge to an incorrect result. When you obtain such an incorrect result, please re-try to deblur with a slightly changed parameters (e.g., using large blur kernel sizes or gamma correction (2.2)). Should you have any questions regarding this code and the corresponding results, please contact Jinshan Pan (sdluran@gmail.com), Zhe Hu (zhu@ucmerced.edu), Zhixun Su (zxsu@dlut.edu.cn) or Ming-Hsuan Yang (mhyang@ucmerced.edu).

