.. File Appendices.rst

Appendices
==========

.. _chap_psfformat:

:file:`.psf` file format description
------------------------------------

PSF models are `FITS binary tables
<https://archive.stsci.edu/fits/fits_standard/node67.html>`_ with a single row,
containing PSF image components. PSF files derived from |MEF|_ images are
themselves |MEF| file, with one ``PSF_DATA`` extension per input image
extension.

A sample of a PSF binary table extension header is given below:

.. code-block:: tcsh

  XTENSION= 'BINTABLE'           / THIS IS A BINARY TABLE (FROM THE LDACTOOLS)    
  BITPIX  =                    8 /                                                
  NAXIS   =                    2 /                                                
  NAXIS1  =                25000 / BYTES PER ROW                                  
  NAXIS2  =                    1 / NUMBER OF ROWS                                 
  PCOUNT  =                    0 / RANDOM PARAMETER COUNT                         
  GCOUNT  =                    1 / GROUP COUNT                                    
  TFIELDS =                    1 / FIELDS PER ROWS                                
  EXTNAME = 'PSF_DATA'           / TABLE NAME                                     
  LOADED  =                 7256 / Number of loaded sources                       
  ACCEPTED=                 6637 / Number of accepted sources                     
  CHI2    =           1.17986252 / Final reduced Chi2                             
  POLNAXIS=                    2 / Number of context parameters                   
  POLGRP1 =                    1 / Polynom group for this context parameter       
  POLNAME1= 'X_IMAGE '           / Name of this context parameter                 
  POLZERO1=   2.046996443272E+03 / Offset value for this context parameter        
  POLSCAL1=   4.074312777519E+03 / Scale value for this context parameter         
  POLGRP2 =                    1 / Polynom group for this context parameter       
  POLNAME2= 'Y_IMAGE '           / Name of this context parameter                 
  POLZERO2=   2.047786695004E+03 / Offset value for this context parameter        
  POLSCAL2=   4.077873875618E+03 / Scale value for this context parameter         
  POLNGRP =                    1 / Number of context groups                       
  POLDEG1 =                    3 / Polynom degree for this context group          
  PSF_FWHM=           2.23724842 / PSF FWHM in image pixels                          
  PSF_SAMP=           0.47601029 / Sampling step of the PSF data                  
  PSFNAXIS=                    3 / Dimensionality of the PSF data                 
  PSFAXIS1=                   25 / Number of element along this axis              
  PSFAXIS2=                   25 / Number of element along this axis              
  PSFAXIS3=                   10 / Number of element along this axis              
  TTYPE1  = 'PSF_MASK'           / Tabulated PSF data                             
  TFORM1  = '6250E   '                                                            
  TDIM1   = '(25, 25, 10)'                                                        
  END                                                                             

The file content is largely self-describing. Please note that

* The ``TDIM1`` keyword in the extension header contains the total number of
  components (conditioned by the ``PSFVAR_DEGREES`` configuration parameter) and
  the dimensions of the tabulated PSF (see the ``PSF_SIZE`` configuration
  parameter).

* Every context parameter (e.g., :math:`X_j`) is rescaled before being used as
  a polynomial variable (:math:`x_j`):

.. math::
  :label: context_rescaling

   x_j = \frac{X_j - \textrm{POLZERO} j}{\textrm{POLSCAL}j}

* The PSF FWHM reported by the ``PSF_FWHM`` keyword is actually derived from
  the mode of the
  `half-light diameter <https://en.wikipedia.org/wiki/Effective_radius>`_
  distribution of the input source image sample.

.. include:: keys.rst

