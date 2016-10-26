.. File Introduction.rst

Introduction
===================

|PSFEx|_ (PSF Extractor) is a computer program that extracts precise models
of the `Point Spread Functions (PSFs)
<http://en.wikipedia.org/wiki/Point_spread_function>`_ from images processed by
|SExtractor|_ and measures the quality of images. The generated PSF
models can be used for model-fitting photometry or morphological
analyses. The main features of |PSFEx| are:

*  Modelling of any arbitrary non-parametric or parametric,
   bandwidth-limited, PSF.
*  Reconstruction of PSF from undersampled images using super-resolution
   on the pixel basis, the Gauss-Laguerre basis :cite:`2005MNRAS.363..197M`
   or a user-provided vector basis.
*  Modelling of PSF variations as a polynomial function of position in
   image, any |SExtractor| measurement, or any numerical |FITS|_ parameter.
*  Tracking of hidden PSF dependencies using `Principal Component
   Analysis <http://en.wikipedia.org/wiki/Principal_component_analysis>`_
   :cite:`2008arXiv0810.0027J`.
*  Computation of PSF homogenisation kernels (to convert variable
   instrumental PSFs to constant round `Moffat
   <http://en.wikipedia.org/wiki/Moffat_distribution>`_
   :cite:`1969A&A.....3..455M` profiles).
*  Automatic selection of point sources.
*  Compatibility with |SExtractor| |FITS| or Multi-Extension |FITS| catalogue
   format in input,
*  |VOTable|_-compliant |XML|_ output of meta-data.
*  |XSLT|_ filter sheet provided for convenient access to metadata from a
   regular web browser.

.. include:: keys.rst

