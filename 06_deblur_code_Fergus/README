README file for Deblurring code used in SIGGRAPH '06 paper: "Removing
Camera Shake from a Single Image" by R.Fergus, B.Singh, A. Hertzmann,
S. T. Roweis, W. T. Freeman.


Distribution code Version 1.2 -- 10/23/2006 by Rob Fergus (fergus at
csail dot mit dot edu) Copyright 2006, Massachusetts Institute of
Technology.

IMPORTANT DISCLAIMER: The source code is messy research-grade code. There
are limited comments through out and lots of strange obsolete sections
and so on as a result of the organic development process. If you are
having problems understanding what the code is doing, please take some
time to work through it yourself and understand what is going
on. Try reading the various references in Section 3 below, hopefully this
should offer some enlightenment. 

The algorithm is not yet sufficiently robust to be a tool that can
used by anyone. You may have to run the algorithm several times on an
image, trying different subwindows etc. before you get a sharp image
out. Experience of using the algorithm certainly does help in
selecting the various input parameters. We are working on improving
the algorithm to bring it up to commerical-grade robustness. Please
remember that blind deconvolution of images is a very challenging
problem that has been around for several decades.

I don't want to be mean but I'd really appreciate it if you could only
email me if you are really stuck. I will try and get back to you, but I may take some time
as I am busy with other stuff.

--------------------------------------------------------------------------------
¦ Contents                                                                     ¦
--------------------------------------------------------------------------------

1. Legal blurb regarding the patent filed on this algorithm.
2. Acknowledgements
3. Academic use 
4. Installing the code
5. Overall structure of algorithm
6. Detailed breakdown of parameter script 
7. Hints and tips for using the algorithm.

-------------------------------------------------------------------------------
¦ 1. Legal blurb                                                              ¦
------------------------------------------------------------------------------- 
Begin patent blurb from MIT TLO

      1. A patent has been filed for this algorithm.
      2. It may be used for non-commercial use only.
      3. If you wish to use it commercially, please contact the MIT
         Technology Licensing Office (tlo@mit.edu): 
      
       Technology Licensing Office
       Massachusetts Institute of Technology
       Five Cambridge Center, Kendall Square
       Room NE25-230
       Cambridge, MA 02142-1493

       Tel: +1 617-253-6966
       Fax: +1 617-258-6790
      
End patent blurb from MIT TLO

-------------------------------------------------------------------------------
¦ 2. Acknowledgements                                                          ¦
------------------------------------------------------------------------------- 

We'd like to acknowledge a variety of people for their help on this
project:
	- David Mackay and James Miskin for graciously allowing us to
	distribute modified versions of their Variational deblurring
	algorithm.
	- Yair Weiss for the Poisson reconstruction code.
	- Antonio Torralba, Don Geman, Marshall Tappen, Fredo Durand
	and others for their useful insights and suggestions.
	- Daniel Dardani and the MIT TLO for their assistance in
	patenting and licensing issues regarding this algorithm.
	- Various funding agencies

-------------------------------------------------------------------------------
¦ 3. Academic use of the code                                                  ¦
------------------------------------------------------------------------------- 

If you use this code to generate results that you use in an academic
publication you should include the following citations in your paper
(included in Bibtex format for your convenience):

- The Fergus et al. SIGGRAPH 06 paper:

@article{Fergus06,
author = "Fergus, R. and Singh, B. and Hertzmann, A. and Roweis, S. T. and Freeman, W.T.",
title =  "Removing Camera Shake From A Single Photograph", 
volume = "25",
issue = "3",
pages = "787--794",
journal = "ACM Transactions on Graphics, {SIGGRAPH} 2006 Conference Proceedings, Boston, {MA}",
publisher = "ACM Press",
year = "2006"
}

- At least one of the following Miskin and Mackay papers:

@incollection{Miskin2000,
author = {James Miskin and David J. C. MacKay},
title = {Ensemble {L}earning for {B}lind {I}mage {S}eparation and {D}econvolution},
booktitle = {Adv.~in {I}ndependent {C}omponent {A}nalysis},
editor = {M. Girolani},
publisher = {Springer-Verlag},
year = {2000}
}

@phdthesis {Miskin2000a,
        author = "Miskin, James W.",
        title = {Ensemble Learning for Independent Component Analysis},
        school = {University of Cambridge},
        year = {2000},
        month = {December},
        url = {http://www.inference.phy.cam.ac.uk/jwm1003/} }

@misc{Miskin_software,
author = {James Miskin},
title = {Train Ensemble Library},
howpublished = "\url{http://www.inference.phy.cam.ac.uk/jwm1003/train_ensemble.tar.gz}",
year = 2000
}



-------------------------------------------------------------------------------
¦ 4. Installing the code                                                       ¦
------------------------------------------------------------------------------- 

4.1 Platform issues / requirements

    - The current implementation is purely in Matlab, so an installation of
      Matlab, along with the Image Processing Toolbox will be
      needed. I've tested the code under Matlab 7.2 and it
      checks out OK., but it should work OK under Matlab 6.5 although
      the included .mat files are V7 format. But you should also read comment 8.6.

    - If you haven't got Matlab then I'm afraid you
      can't use the code. Try getting hold of the Matlab student version.
      
    - It should work under Windows/Mac/Linux, although I've only tried
      it under Linux.

    - As for hardware issues, any modern machine should be fine. The
      memory requirements of the kernel inference part are minimal, although running
      Matlab's deconvlucy on large (>6 megapixel) images can require a
      lot (>1Gb) memory. The inference is slow, so a fast CPU is a
      good idea. I don't think it can be parallelized easily.

    - Hard disk usage: the .mat files generated for each image can be
      quite large due the unnecessary storage of all variables. So you
      might want to have a Gig or two free if you want to run on lots
      of images.
   
4.2 Unpacking and installing the code 

    - The code, prior files and demo images and results are all in one
      .zip file, deblurring.zip (available from
      http://people.csail.mit.edu/fergus)

    - Download this and unzip it somewhere sensible. You should see
      the subdirs: code/ images/ priors/ results/. 
      
    - Ensure that the code/ subdirectory is in your Matlab path.        
  
4.3 Files not included in the .zip file

    - There are certain options in the code (currently disabled)
    that require additional software packages. These are:
    
	- Eero Simoncelli's matlabPyrTools routines
	  http://www.cns.nyu.edu/~lcv/softare.html

-------------------------------------------------------------------------------
¦ 5. Overall structure of the algorithm                                       ¦
------------------------------------------------------------------------------- 

5.1 Directory structure

There are 4 subdirectories. code/ holds all the Matlab source code (it
should be in your Matlab path). priors/ holds the .mat files
specifying the parameters of the mixture of Gaussians priors on image
gradients. images/ holds the original images to be deblurred. Three
example images are provided - add your own in to the same directory.
results/ holds the main input and output files for the algorithm,
detailed below.

5.2 File structure

5.2.1 code/ 

The top-level function is deblur.m This is called with the name of an
image script (passed as a string). This runs the multi-scale inference
procedure to infer the blur kernel. See section 5.3 for
a more detailed explanation.

The core inference algorithm uses source code adapted from James
Miskin and David Mackay who have kindly permitted the redistribution of
the modified versions of their code. All their files start with
'train_'. They are very dense and most likely unintelligible until you
have read the relevant parts of James Miskin's thesis (see Bibtex
reference in section 3.) Please see the header comments in each file
for more details.

fiddle_lucy3.m is the routine that runs Richardson-Lucy and can be set
to save the resulting deblurred output to disk. See comments
at top of file for the various options.

fiddle_lucy4.m - same as fiddle_lucy3.m but runs only on the
intensity channel. For small blurs, this can help reduce the amount
of chrominance noise.

deconvlucy_intens.m runs deconvlucy on intensity channel in YIQ
colorspace. Used by fiddle_lcuy4.m

initialize_parameters2.m is the routine the initializes the gradient
image and kernel estimate at each scale. 

move_level.m does the upsampling after inference at each scale.

GaussianMixtures1D.m, mix_exponentials.m are used to estimate the
prior parameters via EM. They are called by estimate_priors2.m (this
file is a real mess, sorry).

plotgray.m is a useful utility to inspect progress during
inference.

reconsEdge3.m and invDel2.m are two files written by Yair Weiss that do a Poisson
reconstruction from image gradients. Thanks to him for letting me use
it.

The remaining files are just various utilities.

5.2.2 priors/ 

linear_street_4.mat - default priors file. 4 component mixtures of
Gaussians prior estimated from a typical street scene image using
linearized intensities.

linear_whiteboard_4.mat - alternate priors file.  4 component mixtures of
Gaussians prior estimated from an image of a lab whiteboard with
text/equations on it, image using linearized intensities.

5.2.3 images/ 

ian1.jpg - blurry bird picture used in SIGGRAPH 06 presentation

lyndsey2.jpg - blurry statue picture used in SIGGRAPH 06 presentation

5.2.4 results/

ian1.m - image script file for ian1.jpg
ian1_final.jpg - deblurred output from ian1.jpg
ian1_kernel.eps - inferred kernel from ian1.jpg

lyndsey2.m - image script file for lyndsey2.jpg
lyndsey2_final.jpg - deblurred output from lyndsey2.jpg
lyndsey2_kernel.eps - inferred kernel from lyndsey2.jpg

When running the code you will see tmp_<script name>.mat files
appearing in this directory. These may be deleted once inference has
finished.

5.3 How to use the code

Each image you deblur should have a separate image script file (which
is a Matlab .m file). Running this generates a .mat file with the same
name holding all parameters and input variables that the algorithm
will use. The script file is designed to ensure the repeatability of
each experiment. 

Having generated the .mat file, deblur.m is called with the filename
of the .mat file passed in as the only input. This calls the inference
and save the estimate of the kernel and the deblurred image into the
same .mat file. The results may be viewed using the "plotgray"
function.

Here is a run-though of the algorithm from start to end:

1 Copy your blurred image into images/ (e.g. rob1.jpg)
2. Make a copy of one of the example image scripts in results/
   (e.g. 'cp ian1.m rob1.m' if you are using Linux)
3. Edit new image script (e.g. rob1.m), changing the following settings:
	- obs_im to reflect new file name (e.g. obs_im = '../images/rob1.jpg';)
	- AXIS to reflect the area of the image to use for
	inference. You may want to load up in the image in Matlab,
	display it and zoom in to find a suitable region. You can get
	the current axis limits with the command
	'round(axis)'. Cut'n'paste the values into the AXIS setting.
	-  Set BLUR_KERNEL_SIZE by finding  a blurred edge somewhere in the image and estimating
	the size of the blur. Always be generous - its better to give
	a larger value than a smaller one. Minimum sensible values
	would be around 11 or so. Max. sensible values 50-60 or so. 
	- Set initial orientation (FIRST_INIT_MODE_BLUR) by examining
	the rough direction of the blur in the image (horizontal or
	vertical). This one isn't too important. If you're not sure,
	you can specify a delta function.
	- If the image is very large you might want to reduce it by
	setting PRESCALE to be less than 1.
4. Run the image script at your Matlab prompt (e.g. 'rob1'), while in
   the results/ subdirectory. This should generate rob1.mat 
5. Run main inference algorithm with "deblur(<script_name>)", e.g. "deblur('rob1')"
6. Check on progress as it runs by examining the tmp_<script_name>.mat
file using "plotgray(<script_name>)", e.g. "plotgray('rob1')".
7. When it has finished the final results may be viewed with the same
"plotgray(<script_name>)" command. If called with an additional flag
parameter ("plotgray(<script_name>,1)") it will save the kernel and deblurred image out to disk. 
8. If you want to play around with the parameters of the
Richardson-Lucy algorithm, for example using more iterations, or
dropping down a resolution level, then use the "fiddle_lucy3"
function, e.g. "fiddle_lucy3('rob1',20)".

-------------------------------------------------------------------------------
¦ 6. Details of the parameters script		                              ¦
------------------------------------------------------------------------------- 

The script contains all the parameters used by the algorithm. This is
to ensure that each run is entirely repeatable. We now go through all
the different parameters listing what they do.

6.1 User specified parameters

As mentioned, the algorithm isn't quite fully automatic, so you do
need to set some parameters by hand. Typically these are fairly easy
to set once you get the hang of things.

6.1.1 AXIS - 4x1 integer vector. Specifies the subregion of the image
      to use for inferring the kernel. Format is [x_min x_max y_min
      y_max]. These are the coords AFTER any pre-scaling of the image
      with PRESCALE has been applied.
       
6.1.2 BLUR_KERNEL_SIZE - integer.  Size of blur kernel at full
      resolution (after PRESCALE has been applied) in pixels. i.e. 31
      will infer a 31x31 pixel kernel. Although it will probably work
      it is best to use an odd value as it gives a distinct center
      pixels. This parameters indirectly controls the number of scale
      levels to use in inference (along with RESIZE_STEP).

6.1.3 FIRST_INIT_MODE_BLUR - string. Possible values
    {'hbar'¦'vbar'¦'delta'}. Sets the initial 3x3 kernel that will be
    used for inference at the 1st scale level. 'hbar' is a horizontal
    bar (i.e. [0 0 0 ; 1 1 1 ; 0 0 0]/3). 'vbar' is a vertical bar
    (i.e. 'hbar' rotated through 90 degrees). 'delta' is a delta
    function kernel.

6.1.4 obs_im - string. Filename of image to deblur. Can be any format
    that Matlab can handle. If using a RAW file, please convert it
    first to a 16-bit TIFF, and set obs_im to be its filename, and
    then change GAMMA_CORRECTION to 1 (it is assumed that the TIFF has
    values in the range 0 to 65535 since a 1/256 scaling will be
    applied). 

6.1.5 PRESCALE - float in range 0 to 1. Prescaling for
    obs_im. Sometimes it is quicker/easier to work with smaller images
    than the full resolution of obs_im (i.e. if you are 6 megapixels
    or more things can get slow and tedious). By setting this value to
    less than 1, the image will be subsampled using imresize before
    being passed into the algorithm. e.g.  PRESCALE = 0.25 will reduce
    the size of obs_im by a factor of 4. Another time you may want to use
    it is if you have an image with a really huge blur (100+
    pixels). In these cases the inference typically fails, but if you
    subsample the image first, the blur will be of a more manageable
    size.

6.2 Lucy-Richardson parameters

These parameters only relate to the application of the Richardson-Lucy
algorithm to the blurry image, using the inferred kernel. The play no
role in the actual inference stage itself. Since the operation of RL
is somewhat dependent on the shape of the kernel itself, some manual
control may be needed to get the best results. The following
parameters are the most important ones. 

Of course, we are working on doing away with RL altogether, so this
stage of the algorithm should be regarded as a stop-gap measure.

6.2.1 LUCY_ITS - integer. Default = 10. Sensible values 5 to
    30. Number of Lucy-Richardson iterations to use when unblurring
    the image using the inferred kernel. 10 is the default but if the
    kernel is very long and thin you might need more. If you turn it
    up too much it amplifies the noise in the image.

6.2.2 SCALE_OFFSET - integer. Default = 0. Sometimes it may not be
      possible to deblur the image at full resolution due to high
      noise levels (see discussion in 7.3). If this is the case, you
      can tell the fiddle_lucy3 function to use a coarser scale. The
      value of SCALE_OFFSET dictates how many scale levels you drop
      down before selecting the kernel. i.e. SCALE_OFFSET = 1 will use
      a kernel that is one scale level (RESIZE_STEP smaller) than the
      full resolution kernel.

6.2.3 KERNEL_THRESHOLD - float. Default = 7, sensible range is 5 to
    15. This is a threshold on kernel intensities, applied after the
    whole inference stage, designed to remove noise from the
    kernel. It is a dynamic threshold which is the percentage of the
    maximum value of the kernel. Since it is a bit dependent on the
    intensity profile of the kernel, some manual tuning may be needed
    to get the optimal value. If you don't want to use the threshold at all, set it to a
    large number (i.e. 100). If you set it too high, it will start to
    erode the structure of the kernel. This parameter is a bit of a
    hack - in reality you shouldn't need it, but it can make quite a
    bit of difference if the inferred kernel is noisy. 

6.3 Semi-fixed parameters

Occasionally if you are getting no sensible results despite playing
  with the parameters in 6.1 you might want to try changing a few
  other parameters to see if they help. These are listed below in
  rough order of importance.

6.3.1 PRIOR_TYPE - string. Possible values {'street' ¦ 'whiteboard' }.
      Prior file to use. In priors/ are two files with the parameters
      of a mixture of Gaussians fit to the image gradients from
      different scenes. 'street' selects the parameters estimated from
      a typical street scene and should be used by default for most
      images. If you want to use a region of the image which is mainly
      text (i.e. binary valued in nature) then you might want to try
      the 'whiteboard' option. These parameters were estimated from an
      image of a whiteboard with text written on in dark colors. Note
      that both sets of parameters were estimated using *linear*
      gradients (i.e. once gamma correction was removed).

6.3.2 SPECIFIC_COLOR_CHANNEL - integer {0,1,2,3}. Allows you to run
    inference using all 3 or just one of the colors channels. The
    default is 0 which converts the RGB image to grayscale using a
    slightly modified version of rgb2gray (one that understands about
    saturation). 1 selects just the red channel. 2 selects just the
    green channel, 3 selects just the blue channel. This option might
    be useful in night shots when the red and green channels may often
    be saturated but the blue one may not be.

6.3.3 GAMMA_CORRECTION - float. Default = 2.2. Value of gamma used in
    gamma correction. If a RAW file is being deblurred, please change
    this to 1.0.

6.3.4 SATURATION_THRESHOLD - integer 0-255. Intensity value above
    which pixels will be masked out by the saturation mask
    (i.e. not used in inference).

6.3.5 CAMERA_TYPE - string. Default = 'unknown'.  If you have computed
    the tone response map for a specific camera (they are unique to
    each model of camera) then this option allows you to use this
    response map to linearize the intensities of the image. No such
    files are provided with the distribution, but you may wish to
    build your own.

6.3.6 AUTOMATIC_PATCH - binary value. I have written a crude algorithm
    that uses a variety of heuristics to automatically select the
    region of the image to use for inferring the kernel. Setting this
    parameters to 0 turns it off (so the AXIS value above is
    used). Setting it to 1 turns the automatic selection on. Default
    if 0 (off).

6.3.7 SYNTHETIC - binary value. Default = 0. Lets the algorithm know
    that the image provided is actually not blurred and so it should
    synthetically blur it. This was used in testing, but should be
    left set a 0 for all real images. 

6.3.8 BLUR_LOCK - binary value. Default =1. Turns off/on blur kernel
    priors. If off, inference algorithm of Miskin and Mackay will try
    to infer the kernel prior. If on, it uses the parameters
    previously estimated. If you kernel is super-weird then turning this
    *might* help, but its not likely.

6.3.9 CENTER_BLUR - binary value. Default = 1. Turns off/on the
    automatic centering of the blur kernel estimate after inference at
    each scale. This ensures most of the kernel mass is nice and
    central. Occasionally you might need to turn this off if you have
    a long thin tail to the kernel that carries little mass, but I
    recommend just increasing the size of BLUR_KERNEL_SIZE.

 
6.4 Fixed parameters

There are quite a few parameters that should be left fixed to the
values they are set at in the example scripts provided. If I have the
time/energy at some point I will document them individually, but for
now just leave anything below the "Fixed parameters" comment alone.

    
-------------------------------------------------------------------------------
¦ 7. Hints and Tips					                      ¦
------------------------------------------------------------------------------- 

7.1 Saturation
    A bad, bad thing that screws up the algorithm. Try to select a
    region without any saturation (i.e. pixels over 250 in intensity
    or so). There is a mask that is supposed to take out saturated
    pixels but I think there may be some bug with it.  You can
    sometimes avoid saturation by using only one of the color channel
    (see the SPECIFIC_COLOR_CHANNEL parameters).

    Richardson-Lucy also can't handle saturation - it will produce
    nasty color fringes.

7.2 Kernel size
    Always make it quite a bit bigger (30-50% or so) than the guessed
    kernel size. Algorithm doesn't seem too sensitive to the exact
    size as long as it is larger than the kernel.

7.3 Picking the region
    The algorithm is somewhat sensitive unfortunately to the exact
    region picked. Some rough guidelines to picking a good region include:
    - Avoid saturation
    - Areas with nice structured edges
    - Point sources can be used provided they aren't saturated. 
    - Text regions can also be good but change the PRIOR_TYPE to 'whiteboard'
    - Size of the region: larger is usually better but it takes longer.
    - See the examples in the SIGGRAPH paper and slides.

7.4 Initialization
    Algorithm doesn't seem to be too sensitive to this. If you aren't
    too sure just set FIRST_INIT_MODE_BLUR to 'delta'.

7.5 Image noise
   This can be a real problem and is something on which we are
   currently working. This problem wasn't mentioned in the paper as it
   didn't become apparent until we had applied that algorithm to a
   wide range of images.

   If the image is very noisy then the
   noise model tends to breakdown at the fine scale (i.e. towards end
   of the scale sequence) and the kernel can "evaporate" into
   noise. This is due the simple Gaussian model not being able to
   handle the blocky jpeg artifacts. Intuitively, at coarse scales the
   downsampling removes the effects of noise but at fine scales, the
   jpeg block boundaries become visible and generate spurious image
   gradients, i.e. that don't belong to any actual structure within
   the image, confusing the algorithm. HENCE, YOU MAY NOT BE ABLE
   TO DEBLUR THE IMAGE AT FULL RESOLUTION.

   There isn't really any sensible fix for this. One easy option is to
   subsample the image with PRESCALE before starting. You won't have
   such a high-res output, but at least you will have
   something. Alternatively, run using the full resolution and then
   use the SCALE_OFFSET variable to deblur at whatever resolution the
   algorithm worked up until.

   Images from high quality cameras tend to have lower noise than
   cheap devices due to the large sensor area and pixel sizes. Hence
   images from DSLR's usually have no noise problems but those from
   cell phones can be problematic.

   Images taken at night are problematic for noise reasons: the
   automatic grain control on most cameras ramps right up to maximum
   gain and tends to give horrendous noise in the dark areas of the
   image, while the bright areas all tend to be saturated. This
   usually leaves very little of the image which can be used to infer
   the kernel.
   
7.6 Checking on progress
    When you are running the code, it can takes some time. Depending
    on the size of the area selected it could be anywhere from 5
    minutes to an hour or so. I do recommend using the plotgray
    function to examine the progress of the algorithm.

-------------------------------------------------------------------------------
¦ 8. Frequently asked questions
------------------------------------------------------------------------------- 

8.1 Q: Richardson-Lucy Where can I find an explanation of it?
    A: Matlab's deconvlucy is an adaptive damped version of the
    algorithm. An explanation of something pretty close to it can be
    found at: http://www.stsci.edu/stsci/meetings/irw/proceedings/whiter_damped.dir/whiter_damped.html

8.2 Q: How do I generate and run on synthetic images?    
    A: The code used to be able to run in synthetic mode by changing
    SYNTHETIC to 1 but I haven't used this feature in a while and so
    it may need a little debugging. In this mode a sharp image is
    manually blurred and then deblurred again. Although this is a good
    sanity check in that you can successfully unblur images that obey
    your formation model, it doesn't tell you how it will work on real
    images. These contain all kinds of weird non-linearities that can
    screw over your algorithm. Too many papers on blind-deconvolution
    only test on synthetic data. In my opinion, who cares about
    unblurring synthetically blurred images? You want to unblur real
    images, for which there is no sharp version in existence.

8.3 Q: Richardson-Lucy is so old and dumb - it has no image prior! Why
    don't you use something better?
    A: We tried lots of different things. They didn't do significantly
    better than RL and typically took much longer. Many approaches do
    just great on synthetic cases but fall over badly on real data.
    However, we have been working on this and hope to publish
    something soon.

8.4 Q: Why is there so much noise in deblurred image?  
    A: Richardson-Lucy doesn't have an
    explicit image prior (it has an implicit prior by damping the
    reconstruction in an attempt to reduce noise, but it isn't a real
    image prior). Running RL for more iterations makes the image
    sharper but increases the noise levels and the ringing
    artifacts. One particularly noticeable aspect of the noise is the
    colored speckles. These may be removed fairly effectively by
    applying a median filter to the chrominance channels of the
    image. I've not included this in the code since the aim of this
    code distribution is to show the raw output of the algorithm.
   
8.5. Q: What about these ringing artifacts in the deblurred image?
     A: See 8.3 and 8.4. We have been working on it and will publish
     something soon.

8.6  Q: Is the algorithm deterministic?
     A: Yes, given an identical image script file and input image. I
     found that I get very slightly different results if I
     use Matlab 6.5sp1 vs Matlab 7.2. Bizzarrely, in Matlab 6.5, the
     uint8() function acts as floor(), while in Matlab 7.2 it acts as
     round(). Weird.

8.7  Q: You algorithm does not work on my images!
     A: Please understand this raw research: there are many problems
     and limitations of the algorithm as it stands. We are currently
     working to make it more robust and extend its applicability. 
     But is important to understand the limitations of the algorithm,
     at present:
	   1. It can only remove specific types of blur - it won't
     work with out-of-focus blur, or when the camera movement had
     large in plane rotation component (check the blur pattern in two
     opposite corners - if they look very different, you have a large
     in-plane rotation). If you wave the camera around at night with
     the shutter open for 10 seconds, its not going to work either! 
     The practical maximum blurs I've managed to remove are in the 
     range 50-60 pixels across.
	  2. If you have high image noise then it isn't going to
     unblur at the full resolution, see 7.5.
	  3. The algorithm it not totally automatic - most notably,
     it is senstive to the region selected. Try finding regions where
     you can see the blur fairly clearly. Try using a window of at least
     200x200 pixels or so.
	  4. Night images are particuarly hard (see 7.5 again).
	  
8.8  Q: I've written my own deblurring algorithm and it works much better than yours.
     A: Great, you should submit a paper explaning your algorithm and
     showing results on real images to an academic conference. That
     way you can tell everyone about your approach and help the
     community make progress on this challenging problem. If you file
     a patent before you publish your paper, your invention will be protected.


