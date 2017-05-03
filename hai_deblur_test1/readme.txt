This package contains an implementation of the motion blur kernel estimation algorithm described in the paper: Jinshan Pan, Risheng Liu, Zhixun Su, and Xianfeng Gu, "Kernel Estimation from Salient Structure for Robust Motion Deblurring", Signal Processing: Image Communication, 28(9): 1156-1170 (2013)


Run "demo.m" to try the examples included in this package.

Note: 
1. Please note that the algorithm sometimes may converge to an incorrect result. When you obtain such an incorrect result, please re-try to deblur with a slightly changed parameters (e.g., using large blur kernel sizes).

2. Most blurred images in the "image" folder are obtained from the package which is available at http://appsrv.cse.cuhk.edu.hk/~xuli/robust_deblur_executable.zip  


3. If the kernel is not at the center, you can use "adjust_psf_center.m" to adjust it.

The code is tested in MATLAB 2011b(32bit) under the MS Windows 7 64bit version with an Intel Xeon CPU@2.53GHz and 12GB RAM.

Should you have any questions regarding this code and the corresponding results, please contact Jinshan Pan via email (jspan@mail.dlut.edu.cn).

