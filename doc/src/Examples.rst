.. File Examples.rst

Examples
========

In the following, examples of use of PSFEx are given, together with
commented command lines.

Hands-on example 1
------------------

Let us consider a V band FITS image RX\_J2202-19\_V.fits and its weight
map RX\_J2202-19\_V.weight.fits. We wish to fit all the galaxies of the
image with a galaxy model using SExtractor, which requires computing a
model of the PSF first.

We must first run SExtractor on this image to obtain a temporary
catalogue in FITS\_LDAC format that contains small sub-images from which
the PSF model will be extracted. For this, we define in a SExtractor
parameter file — let us call it prepsfex.param — the parameters required
for the use of PSFEx::

    X_IMAGE
    Y_IMAGE
    FLUX_RADIUS
    FLUX_APER(1)
    FLUXERR_APER(1)
    ELONGATION
    FLAGS
    SNR_WIN
    VIGNET(35,35)

TBW

Example 2: very wide photographic plate
---------------------------------------

TBW

Example 3: unfocused instrument
-------------------------------

TBW

Frequently Asked Questions
===========================

   
Skeptical Sam doesn’t have time to test software extensively but is
always keen on asking aggressive questions to the author to find out if
a program could fit his needs.

**PSFEx represents PSFs as an array of tabulated values! Can
it really deal with undersampled images? Isn’t it too noisy?**

PSFEx was designed from the ground up to deal with
undersampled images and arbitrary PSFs. Although the PSF “model” in
PSFEx is actually a small image, it is sampled at a different step than
the original pixels: more finely for undersampled observations, and more
coarsely for oversampled observations, to avoid any loss and redundancy
of information. Despite built-in regularisation, PSF models
reconstructed on the pixel basis can indeed be noisy if the number of
selected stars is small. This can be circumvented to some extent by
using *ad hoc* basis to solve for the PSF model coefficients.

**I heard that PSFEx has been developed almost 12 years ago,
and has been used for production at TERAPIX for many years. Why have you
waited until 2010 for releasing it to the general community?**

PSFEx was originally developed for doing PSF-fitting
crowded-field photometry with SExtractor. However I was not very happy
with the way it worked, as SExtractor’s detection and deblending engine
is not meant to deal with crowded star fields. The current release of
PSFEx is made in the framework of the EFIGI [2]_ and DES [3]_ projects,
as a support tool for galaxy model-fitting.

**I would like to use PSFEx to generate PSF models for
weak-lensing analyses. Is it the right tool for that?**

Simulations of 1h exposures with a 4m optical
telescope and sub-arcsecond seeing show that ellipticities of galaxies with
a Signal-to-Noise Ratio SNR\ :math:`>20` can be recovered with a level of
systematics below :math:`10^{-3}` using PSFEx models, even in the presence
of significant amounts of coma and astigmatism. This is for constant
PSFs. Tests with variable PSFs are ongoing.


Troubleshooting
===============

TBW

