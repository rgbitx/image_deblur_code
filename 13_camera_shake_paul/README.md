Correcting Camera Shake by Incremental Sparse Approximation
=================================================

This package contains a MATLAB implementation of the method described in [1] and results from this method. The code was developed by modifying the L1/L2 blind deconvolution method [2] and employs the SPG solver `minConf_SPG`, available at [3].

(c) Paul Shearer, 2012. Permission to use, copy, and modify this code is given for academic research purposes only. Commercial use without permission is prohibited.

Startup Notes
-------------

- To add the package's subfolders to your path, use the startup_cam script. 
- To avoid clashes between this package's functions and existing functions on your MATLAB path, it may be helpful to remove your working directory from the path before starting.

Key Scripts
-----------

- `./test_blind_deconv`: deconvolves an image provided as a matlab array or as a file (e.g. mukta.jpg), or an image from the Levin test set.
- `./decon_levin_testset`: reproduces the results of [1] on Levin et al's test set of 32 test images.
- `./show_levin_results`: shows the MSE results of our run on Levin's data
- `./benchmark_results_slideshow`: slideshow comparing the three methods in the paper on all 32 test images from the Levin set.

Code flow for `test_blind_deconv`
--------------------------------

    [kernel estimation]
    >> ms_blind_deconv
       >> ss_blind_deconv
         - update x using projected gradients
         - update k using minConf_SPG
       - upsample PSF and image and repeat ss_blind_deconv
    [nonblind deconvolution with estimated kernel]

Data (in `./data`):
-----------------
    im*ker* : Levin's test set (4 images, 8 kernels)
    isep*   : results on Levin test set, our algorithm
    diagfe* : results on Levin test set, Levin et al's marginalization algorithm
    babacan_results.mat : results on Levin test set, Babacan et al's variational Bayes algorithm 


References
----------
[1] Paul Shearer, Anna C. Gilbert, Alfred O. Hero III. Correcting Camera Shake by Incremental Sparse Approximation. [arXiv link, 2013]

[2] D. Krishnan, T. Tay, and R. Fergus, “Blind deconvolution using a normalized sparsity measure,” in IEEE Conference on Computer Vision and Pattern Recognition, 2011, pp. 233–240.

[3] Mark Schmidt, http://www.di.ens.fr/~mschmidt/Software/minConf.html