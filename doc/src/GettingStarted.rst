.. File GettingStarted.rst

.. _gettingStarted:

Getting started
===============

|PSFEx| is run from the shell with the following syntax:

.. code-block:: console

 $ psfex Catalog1 [Catalog2 ...] -c configuration-file [-Parameter1 Value1 -Parameter2 Value2 ...]

The parts enclosed within brackets are optional. The file names of input
catalogues can be directly provided in the command line, or in lists that are
ASCII files with each catalogue name preceded by ``@`` (one per line). One
should use lists instead of the catalogue file names if the number of input
catalogues is too large to be handled directly by the shell. Any
`-Parameter Value` statement in the command-line overrides the 
corresponding definition in the configuration file or any default value (see
below).

Input files
-----------

Catalogues
~~~~~~~~~~

|PSFEx| does not work directly on images. Instead, it operates on |SExtractor|_
catalogues that have a small sub-image ("vignette") recorded for each
detection. This makes things much easier for |PSFEx| as it does not have to
handle the detection and deblending processes. The catalogue files read by
|PSFEx| must be in |SExtractor| FITS_LDAC binary format. This allows
|PSFEx| to have access to the original image header content. The catalogues
*must* contain all the following parameters in order to be processable by
|PSFEx|:

* small image ("vignette") centered on the object, :samp:`VIGNET({w},{h})`,
  where :samp:`{w}` and :samp:`{h}` are respectively the width and the height of 
  the image in pixels,

* centroid coordinates, e.g. ``X_IMAGE`` and ``Y_IMAGE``,

* half-light radius ``FLUX_RADIUS``,
* Signal-to-noise ratio in a Gaussian window ``SNR_WIN``

* flux measured through a fixed aperture, e.g. ``FLUX_APER(1)``,

* flux uncertainty, e.g. ``FLUXERR_APER(1)``,

* object elongation ``ELONGATION``,

* extraction flags ``FLAGS``.

The ``VIGNET`` dimensions :math:`w` and :math:`h` set the maximum size, in
pixels, of the image area stored for each detection. It is advised to
use square sub-images (:math:`w = h`) with an odd number of pixels (for
symmetry across the central pixel), and so that the sub-image covers
much of the visible footprint of non-saturated stars.

The size of sub-images in the catalogue (here :math:`35\times 35` pixels)
must been chosen so that the frame encloses much of the visible
footprint of non-saturated stars. It is recommended not to use
excessively large sub-images as they lead to unpractically large
catalogues and make the |PSFEx| PSF cleaning procedure less robust. In
practice, values such as in ``VIGNET(45,45)`` for instance, will generally
work well with most images.

|SExtractor| configuration settings for the pre-|PSFEx| run do not require
much tuning in general. The |SExtractor| configuration file only need to differ
from the default one on a few keywords. However these keywords must be set with
care:

* ``CATALOG_TYPE`` should be set to ``FITS_LDAC`` (binary |FITS|_ catalogue).

* ``PHOT_APERTURE`` defines the diameter of the circular aperture (in pixels)
  used as a reference for normalising the amplitude of the PSF model. It should
  be set to a value large enough so that variations due to seeing or
  aberrations are negligible at the level required for photometric
  analyses. But be aware that using excessively large apertures lead to noisy
  measurements and are more prone to light pollution by neighbouring
  sources. For professional, ground-based images, a value corresponding to an
  aperture diameter of :math:`5''` is often a good compromise.

* Detector gain: The effective “gain” (or more exactly the conversion factor,
  in units of electrons/ADU) is required by |SExtractor| to compute the
  uncertainty on pixel values, especially with bright star images.
  The ``GAIN_KEY`` configuration parameter tells SExtractor what keyword in the
  original |FITS| image header carries the detector gain. The default string for
  ``GAIN_KEY`` is ``GAIN``. In multi-CCD cameras, the gain can slightly vary
  from one CCD to another, and can also vary with time. When working with
  single exposures, it is recommended to let |SExtractor| read the gain value
  from the image header, if it is present. In some detector chips with multiple
  readouts, several values of the gain may be present in the same header under
  different keywords (e.g. ``GAINA``, ``GAINB``). Since the gain differences
  will often be negligibly small for |PSFEx| (10% or less), it is usually safe
  to use one single value for the whole chip (for example, ``GAINA``).
  Note that what matters for reduced images is the “effective” gain, not the
  original detector gain. Dividing pixel values by some amount (e.g., the
  exposure time or a non-normalised flat field) multiplies the effective gain
  by the same amount. Combining several images also modifies the gain. For
  example if the final image is the mean of several images, the effective gain
  will be equal to the initial gain multiplied by the number of images. Taking
  the median or a weighted average also affect the gain (see, e.g., the
  `SWarp <https://www.astromatic.net/pubsvn/software/swarp/trunk/doc/swarp.pdf>`_
  documentation). Recent versions of |SWarp| properly take into account the
  effect of flux scaling and stacking on the gain, and insert the average
  effective gain in the output image header.
  Finally if no keyword with the name specified by ``GAIN_KEY`` can be found in
  the FITS image header, |SExtractor| will fall back to using the gain value
  specified by the ``GAIN`` configuration parameter. The default fallback value
  is ``0``, which actually tells SExtractor that the gain should be considered
  infinite (bright pixels not noisier that faint pixels).

* Saturation level: |PSFEx| requires |SExtractor| to flag all saturated sources,
  which may otherwise contaminate the “clean” star sample used to compute the
  PSF model. |SExtractor| identifies saturated sources by checking if the value
  of at least one source pixel exceeds a given “saturation level”.
  As for the gain above, |SExtractor| examines first the value of a FITS image
  header keyword to read the saturation level (in ADUs). The header keyword can
  be set with the ``SATUR_KEY`` |SExtractor| configuration parameter; the default
  string for ``SATUR_KEY`` is ``SATURATE``. ``SATURATE`` is commonly found in
  image files released by observatories and pipelines.  Unfortunately, in
  practice, it is often found to be set at a value higher than that at which
  the detector markedly starts to behave non-linearly. It is therefore highly
  recommended to examine *visually* saturated stars on images and check if
  pixel values systematically exceed the saturation level reported in the
  header. If not, it is advised to give ``SATUR_KEY`` a name unlikely to exist
  as a keyword in the FITS file (for example, ``DUMMY``), and force the saturation
  level in |SExtractor| to a lower value using the ``SATUR_LEVEL`` fallback
  parameter.  Note that some detector/amplifier combinations start becoming
  non-linear at levels below the apparent saturation limit, so it is always
  safer to give a saturation level about 10% lower than the lowest value
  derived from the visual examination of all images.

Output files
------------

PSF model files (:file:`.psf`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main purpose of |PSFEx| is to create a PSF model for each of the images from
which the input catalogues were extracted. The PSF models are stored under file
names that are given the :file:`.psf` extension by default (this may be changed
with the ``PSF_SUFFIX`` configuration parameter). The :file:`.psf` files are
FITS binary tables that can be read back into |SExtractor| to perform accurate
model-fitting of the sources being detected. A detailed description of the
:file:`.psf` file format is given :ref:`in the Appendix <chap_psfformat>`.

PSF homogenisation files (:file:`.homo`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is presently an experimental feature. In addition to computing PSF models,
|PSFEx| has the possibility to derive “PSF homogenisation kernels” for all input
catalogues. A PSF homogenisation kernel is a (variable) convolution kernel
which, when applied to an image, gives the point sources it contains a
constant, arbitrary shape. For practical purposes the target shape will
preferably be a perfectly round analytical function, such as a Moffat
:cite:`1969A&A.....3..455M` profile:

.. math::
  :label: moffat
   I(r) = I_0 \,\left [1+\left ({ r\over a} \right )^2 \right ]^{-\beta}

Homogenising the PSF of a set of images can allow for more consistent
image combinations and measurements, once the consequences on noise have
been properly taken into account.

PSFEx stores PSF homogenisation kernels as FITS data cubes. File names are
given the :file:`.homo` extension by default; this may be changed using the
``HOMOKERNEL_SUFFIX`` configuration parameter. :file:`.homo` files can be read
by the :program:`PSFnormalize` software developed by Tony Darnell from the
Dark Energy Survey data-management team to perform fast convolution of the
original images :cite:`2009ASPC..411...18D`. The |SWarp| software may also later
include this possibility.

Diagnostic files
~~~~~~~~~~~~~~~~

Three types of files can be generated by |PSFEx|, providing diagnostics
about the derived PSF and the modelling process:

* "Check-images" are basic FITS files containing images of the PSF model, fit
  residuals, etc.. Configuration parameters ``CHECKIMAGE_TYPE`` and
  ``CHECKIMAGE_NAME`` allow the user to provide a list of check-image types and
  file names, respectively, to be produced by |PSFEx|. A complete list of
  available check-image types is given in §[chap:paramlist]. Many check-images
  are actually aggregates of several small images; they may be stored as grids
  (the default) or as datacubes if the ``CHECKIMAGE_DATACUBE`` parameter is set
  to ``Y``.
* "Check-plots" are graphic charts generated by PSFEx, showing maps or trends
  of PSF measurements. The ``CHECKPLOT_TYPE`` and ``CHECKPLOT_NAME`` configuration
  parameters allow the user to provide a list of check-plot types and file
  names, respectively. A variety of raster and vector file formats, from JPEG
  to Postscript, can be set with ``CHECKPLOT_DEV`` (the default format is
  PNG). See the :ref:`CHECKPLOT <param_checkplot>` section of the
  :ref:`configuration parameter list <param_list>` below for details.
* An |XML|_ file providing a processing summary and various statistics in
  |VOTable|_ format is written if the ``WRITE_XML`` switch is set to ``Y``
  (the default). The ``XML_NAME`` parameter can be used to change the default
  file name :file:`psfex.xml`. The |XML| file can be displayed with any recent
  web browser; the |XSLT| stylesheet installed together with |PSFEx| will
  automatically translate it into a dynamic, user-friendly web-page
  (:numref:`fig_psfexxml`). For more advanced usages (e.g., access from a
  remote web server), alternative |XSLT| translation URLs may be specified
  using the ``XSL_URL`` configuration parameter.

.. _fig_psfexxml:

.. figure:: figures/psfex_xml.*
   :figwidth: 100%
   :align: center

   Rendition of a psfex.xml |XML|-|VOTable| file generated by 
   |PSFEx| with the `Firefox <http://www.mozilla.org/firefox>`_ web-browser.


The Configuration file
----------------------

Each time it is run, |PSFEx| looks for a configuration file. If no
configuration file is specified in the command-line, it is assumed to be
called :file:`default.psfex` and to reside in the current directory. If no
configuration file is found, |PSFEx| will use its own internal default
configuration.

Creating a configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|PSFEx| can generate an ASCII dump of its internal default configuration, using
the ``-d`` option. By redirecting the standard output of |PSFEx| to a file, one
creates a configuration file that can easily be modified afterwards:

.. code-block:: console

  $ psfex -d > default.psfex

and a more extensive dump with less commonly used parameters can be generated
by using the ``-dd`` option.


Format of the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The format is ASCII. There must be only one parameter set per line,
following the form::

 Config-parameter     Value(s)

Extra spaces or linefeeds are ignored. Comments must begin with a ``#``
and end with a linefeed. Values can be of different types: strings (can
be enclosed between double quotes), floats, integers, keywords or
Boolean (`Y/y` or `N/n`). Some parameters accept zero or several values,
which must then be separated by commas. Integers can be given as
decimals, in octal form (preceded by digit 0), or in hexadecimal
(preceded by `0x`). The hexadecimal format is particularly convenient for
writing multiplexed bit values such as binary masks. Environment
variables, written as ``$HOME`` or ``${HOME}`` are expanded.

.. _param_list:

Configuration parameter list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a list of all the parameters known to |PSFEx|. Please refer to
next section for a detailed description of their meaning. Some
“advanced” parameters (indicated with an asterisk) are also listed. They
must be used with caution, and may be rescoped or removed without notice
in future versions.

:Parameter: ``BADPIXEL_FILTER``\*
:Default: ``N``
:Type: *Boolean*

  If true (``Y``), input objects with vignettes containing more than
  ``BADPIXEL_NMAX``  pixels flagged by |SExtractor| as bad or from deblended
  neighbours will be rejected. 

-----------------

:Parameter: ``BADPIXEL_NMAX``
:Default: ``N``
:Type: *Boolean*

  Maximum number of bad pixels tolerated in the vignette before an object
  is rejected (``BADPIXEL_FILTER`` must be set to ``Y``)

-----------------

:Parameter: ``BASIS_NAME``
:Default: ``basis.fits``
:Type: *String*

  File name for the user-supplied FITS datacube of basis vector
  images (``BASIS_TYPE``} must have been set to ``FILE``)


----------------

:Parameter: ``BASIS_NUMBER``
:Default: ``20``
:Type: *integer*

  Size of basis vector set: square-root of the number of pixels for
  ``BASIS_TYPE PIXEL``, :math:`n_{\rm max}` for ``BASIS_TYPE GAUSS-LAGUERRE``, or
  number of vectors for ``BASIS_TYPE FILE``.

--------------

:Parameter: ``BASIS_SCALE`` 
:Default: ``1.0``
:Type: *float*

  Scale size of ``BASIS_TYPE GAUSS-LAGUERRE`` vector images.

-----------------

:Parameter: ``BASIS_TYPE``
:Default: ``BASIS_AUTO``
:Type: *keyword*

  Basis vector set:

  - ``NONE``: No basis; the PSF is derived solely from the robust
  - ``PIXEL``: Pixel basis (super-resolution)
  - ``PIXEL_AUTO``: Equivalent to ``NONE`` for properly sampled images;
    switches automatically to ``PIXEL`` (super-resolution) for critically sampled
    and undersampled data.
  - ``GAUSS_LAGUERRE``: Gauss-Laguerre basis (also known as *polar shapelets* in the weak-lensing community).
  - ``FILE``: User-supplied vector basis, in the form of a FITS datacube (see ``BASIS_NAME``).

----------------

:Parameter: ``CENTER_KEYS``
:Default: ``X_IMAGE, Y_IMAGE``
:Type: *strings*

  Catalogue ``Keys`` (|SExtractor| measurement parameters) that
  define the initial guess for the source coordinates. Note that all input
  ``vignettes`` are automatically re-centred by |PSFEx| using an iterative
  Gaussian-weighted algorithm, hence the centring parameter is not critical.

-------------------------

:Parameter: ``CHECKIMAGE_CUBE`` 
:Default: ``N``
:Type: *Boolean*

  If true (``Y``), check-images will be saved as data-cubes.

-------------------------

:Parameter: ``CHECKIMAGE_NAME``
:Default: ``chi.fits,proto.fits,samp.fits,resi.fits,snap.fits``
:Type: *strings*

  File name of the check-image (diagnostic FITS image) of each type
  (:program:`.fits extension` is not required, as it is assumed by default).

-----------

:Parameter: ``CHECKIMAGE_TYPE``
:Default: ``CHI, PROTOTYPES, SAMPLES, RESIDUALS, SNAPSHOTS``
:Type: *keywords*

  Types of check-images (diagnostic FITS images) to generate during |PSFEx| processing:

   - ``NONE`` No check-image.
   - ``CHI`` (square-root of) :math:`\chi^2` maps for all input vignettes.
   - ``PROTOTYPES`` Versions of input vignettes, recentred, rescaled and resampled to PSF resolution.
   - ``SAMPLES`` Input vignettes in their original position, resolution and flux scaling.
   - ``RESIDUALS`` Input vignettes with best-fitting local PSF models
     subtracted.
   - ``SNAPSHOTS`` Grid of PSF model snapshots reconstructed at each
     position/context.
   - ``MOFFAT`` Grid of :ref:`Moffat models <eq_moffat>` fitted to PSF
     model snapshots at each position/context.
   - ``-MOFFAT`` Grid of PSF model snapshots reconstructed at each
     position/context with best-fitting :ref:`Moffat models <eq_moffat>` 
     subtracted.
   - ``-SYMMETRICAL`` Grid of PSF model snapshots reconstructed at each
     position/context with symmetrised image subtracted.
   - ``BASIS`` Basis vector images used by |PSFEx| to model the PSF.

-------------------------

.. _param_checkplot:

:Parameter: ``CHECKPLOT_ANTIALIAS`` 
:Default: ``Y``
:Type: *Boolean*

  If true (``Y``), PBM, PNG and JPEG check-plots are generated with
  anti-aliasing. `ImageMagick <http://www.imagemagick.org>`_ 's ``convert`` tool must be installed.

-------------------------

:Parameter: ``CHECKPLOT_DEV``
:Default: ``PNG``
:Type: *keywords*

  PLPlot devices to be used for check-plots (all devices may not be available, see PLPlot documentation for details): 

   - ``NULL`` No output
   - ``XWIN`` X-Window
   - ``TK`` Tk window (if available)
   - ``XTERM`` XTerm window
   - ``AQUATERM`` AquaTerm window (Mac OS X)
   - ``PLMETA`` PLPlot .plm} meta-file
   - ``XFIG`` XFig .fig} vector file
   - ``LJIIP`` HP LaserJet IIP .lj} bitmap file
   - ``LJ_HPGL`` HP LaserJet .hpg HPGL vector file
   - ``IMP`` Impress .imp} file
   - ``PBM`` Portable BitMap .pbm image
   - ``PNG`` Portable Network Graphics .png image
   - ``JPEG`` JPEG .jpg image
   - ``PDF`` Portable Document Format .pdf file
   - ``PS`` Black-and-white .ps Postscript file
   - ``PSC`` Colour .ps Postscript file
   - ``PSTEX`` PSTeX (a variant of Postscript) .ps file

-------------------------

:Parameter: ``CHECKPLOT_NAME``
:Default: ``fwhm, ellipticity, counts, countfrac, chi, resi``
:Type: *strings*

  File names for each series of check-plots. |PSFEx| will 
  automatically insert the associated catalogue names, and append/replace
  file name extensions with the appropriate ones,  depending on the chosen
  ``CHECKPLOT_DEV(s)`` (.png for PNG files, .jpg for JPEG, etc.).

-------------------------

:Parameter: ``CHECKPLOT_RES`` 
:Default: ``0`` 
:Type: *integers*  (:math:`n \le 2`)

  Check-plot x,y resolution for bitmap devices (``0`` is equivalent to ``800,600``).

-------------------------

:Parameter: ``CHECKPLOT_TYPE``
:Default: ``FWHM, ELLIPTICITY, COUNTS, COUNT_FRACTION, CHI2, RESIDUALS``
:Type: *keywords*

  Diagnostic check-plots to be generated during |PSFEx| processing (|PSFEx| must  have been configured without the ``--without-plplot option``):    

   - ``NONE`` No plot.
   - ``FWHM`` Map of the model PSF Full-Width at Half-Maximum over the field
     of view (one for each input catalogue). 
   - ``ELLIPTICITY`` Map of the model PSF ellipticity over the field of view
     (one for each input catalogue). 
   - ``COUNTS`` Map of the spatial density of point sources (initially)
     selected over the field of view (one for each catalogue). 
   - ``COUNT_FRACTION`` Map of the fraction of point sources accepted over the
     field of view (one for each catalogue). 
   - ``CHI2`` Map of the average :math:`\chi^2/{\rm d.o.f.}` over the field of
     view (one for each catalogue). 
   - ``MOFFAT_RESIDUALS`` Map of Moffat (eq.~[\ref{eq:moffat}]) residual
     indices over the field of  view (one for each catalogue). 
   - ``ASYMMETRY`` Map of asymmetry indices over the field of view (one for
     each catalogue). 

-------------------------

:Parameter: ``HOMOBASIS_NUMBER`` 
:Default: ``10``
:Type: *integer*

  Size of the homogenisation kernel basis vector set: :math:`n_{\rm max}` for
  ``HOMOBASIS_TYPE GAUSS-LAGUERRE``.

-------------------------

:Parameter: ``HOMOBASIS_SCALE`` 
:Default: ``1.0``
:Type: *float*

  Scale size of ``HOMOBASIS_TYPE GAUSS-LAGUERRE`` homogenisation
  kernel vector images.

-------------------------

:Parameter: ``HOMOBASIS_TYPE`` 
:Default: ``NONE``
:Type: *keyword*

  Basis vector set for the homogenisation kernel:

   - ``NONE`` No basis; no homogenisation kernel is computed.}
   - ``GAUSS_LAGUERRE`` Gauss-Laguerre basis (also known as `polar shapelets`
     in the weak-lensing community). 

-------------------------

:Parameter: ``HOMOKERNEL_SUFFIX`` 
:Default: ``.homo.fits``
:Type: *string*

  Filename suffix of the homogenisation kernels computed by |PSFEx|.

-------------------------

:Parameter: ``HOMOPSF_PARAMS}``\*
:Default: ``2.0, 3.0``
:Type: *floats* (:math:`n \le 2`)

  :ref:`Moffat model <eq_moffat>` Full-Width at Half-Maximum and :math:`\beta` parameters
  of the idealised target PSF chosen for homogenisation.

-------------------------

:Parameter: ``MEF_TYPE`` 
:Default: ``INDEPENDENT``
:Type: *keyword*

  How |PSFEx| should deal with multi-extension catalogues (extracted from mosaic camera images):
   - ``INDEPENDENT`` Derive the PSF model for each extension independently.
   - ``COMMON`` Derive a common PSF model for all extensions.

-------------------------

:Parameter: ``NEWBASIS_NUMBER`` 
:Default: ``8``
:Type: *integer*

  Size of the image vector set (number of basis vectors) derived  by
  |PSFEx| from the input vignettes.}

-------------------------

:Parameter: ``NEWBASIS_TYPE`` 
:Default: ``NONE``
:Type: *keyword*

  Type of image vector bases derived from input vignettes by |
   - ``NONE`` No basis is computed.
   - ``PCA_MULTI`` Karhunen-Lo\`eve basis from Principal Component
     Analysis on all FITS extensions.
   - ``PCA_SINGLE`` Karhunen-Lo\`eve bases from Principal Component
     Analysis on individual FITS extensions.

-------------------------

:Parameter: ``NTHREADS``
:Default: ``0``
:Type: *integer*

  Number of threads (processes) to be used for parallel computation.
  |PSFEx| must have been configured with the ``--disable-threads`` option at
  compile time for this parameter to take effect. Note that multi-threading is
  disabled in the current version of |PSFEx|

-------------------------

:Parameter: ``PHOTFLUX_KEY``
:Default: ``FLUX_APER(1)``
:Type: *string*

  Catalogue ``Key`` (|SExtractor| measurement parameter) that
  defines the flux of sources, and therefore the normalisation of the PSF
  amplitude. It is recommended to use a fixed aperture magnitude; the aperture
  diameter set in |SExtractor| should be large enough so that the fraction of
  flux enclosed stays constant from point source to point source, and small enough
  to preserve the signal-to-noise ratio.

-------------------------

:Parameter: ``PHOTFLUXERR_KEY``
:Default: ``FLUXERR_APER(1)``
:Type: *string*

  Catalogue ``Key`` (|SExtractor| measurement parameter) that
  defines the flux measurement uncertainty on each source. It is used for
  computing the source signal-to-noise ratio.

-------------------------

:Parameter: ``PSF_ACCURACY``
:Default: ``0.01``
:Type: *float*

  Expected accuracy of vignette pixel values (standard deviation of
  the flux fraction).

-------------------------

:Parameter: ``PSF_PIXELSIZE``
:Default: ``1.0``
:Type: *float*

  Effective pixel size (width of the top-hat intra-pixel response
  function) in pixel step units.

-------------------------

:Parameter: ``PSF_RECENTER``
:Default: ``Y``
:Type: *Boolean*

  If true (``Y``), input vignettes are recentred at each iteration
  of the PSF modelling process.

-------------------------

:Parameter: ``PSF_SAMPLING``
:Default: ``0.0``
:Type: *float*

  Sampling step of the PSF models, in pixels. Use ``0`` for automatic
  sampling.

-------------------------

:Parameter: ``PSF_SIZE``
:Default: ``25, 25``
:Type: *integers* (:math:`n \le 2`)

  Dimensions of the tabulated PSF models, in PSF ``pixels``.

-------------------------

:Parameter: ``PSF_SUFFIX`` 
:Default: `` .psf``
:Type: *string*

  Filename suffix for PSF models computed by |PSFEx|.

-------------------------

:Parameter: ``PSFVAR_DEGREES``
:Default: ``2``
:Type: *integers* (:math:`n = n_{\rm groups}`)

  Degree of  polynomial of each context group. ``0`` indicates a
   constant PSF.

-------------------------

:Parameter: ``PSFVAR_GROUPS``
:Default: ``1, 1``
:Type: *integers* (:math:`n = n_{\tt PSFVAR\_KEYS}`)

  Polynomial group which each context key belongs to.

-------------------------

:Parameter: ``PSFVAR_KEYS``
:Default: ``X_IMAGE, Y_IMAGE``
:Type: *strings* (:math:`n \le 2`)

  List of ``keys`` (|SExtractor| measurement parameters) on which
  the PSF is supposed to depend (e.g. ``X_IMAGE, Y_IMAGE`` for a spatial
  mapping of the PSF). Keywords preceded with a colon are interpreted as FITS
  image keywords instead of |SExtractor| parameters.

-------------------------

:Parameter: ``PSFVAR_NSNAP`` 
:Default: ``9``
:Type: *integer*

  Number of PSF snapshots computed on each axis. This also defines the
  resolution of the grid on which diagnostics and check-plot maps are computed.

-------------------------

:Parameter: ``SAMPLE_AUTOSELECT``
:Default: ``Y``
:Type: *Boolean*

  If true (``Y``), input vignettes are automatically selected based
  on the source FWHMs, inside the range specified by ``SAMPLE_FWHMRANGE``,
  with fractional FWHM variability ``SAMPLE_VARIABILITY``.

-------------------------

:Parameter: ``SAMPLE_FLAGMASK`` 
:Default: ``0x00fe``
:Type: *integer*

  Bit mask applied to SExtractor flags for rejecting input vignettes.

-------------------------

:Parameter: ``SAMPLE_FWHMRANGE``\*
:Default: ``2.0, 10.0``
:Type: *floats* (:math:`n=2`)

  Range (in pixels) of source FWHMs (Full-Width at Half-Maximum) allowed
  for input vignettes. FWHMs are currently estimated based on |SExtractor|'s
  ``FLUX_RADIUS`` measurements.

-------------------------

:Parameter: ``SAMPLE_MAXELLIP``
:Default: ``0.3``
:Type: *float*

  Maximum source ellipticity (i.e. :math:`\frac{\tt{A\_IMAGE} - \tt{B\_IMAGE}}{\tt{A\_IMAGE} + \tt{B\_IMAGE}}`) allowed for input vignettes.

-------------------------

:Parameter: ``SAMPLE_MINSN``
:Default: ``20.0``
:Type: *float*

   Minimum source Signal-to-Noise ratio allowed for input vignettes.

-------------------------

:Parameter: ``SAMPLE_VARIABILITY``
:Default: ``0.2``
:Type: *float*

  Maximum fractional FWHM variability (1.0 = 100%) allowed for input vignettes.

-------------------------

:Parameter: ``SAMPLEVAR_TYPE`` 
:Default: ``SEEING``
:Type: *keyword*

  Catalogue-to-catalogue variability criteria for vignette selection:

  - ``NONE`` No differences between catalogues.
  - ``SEEING`` Seeing (hence FWHM) is expected to vary.

-------------------------

:Parameter: ``STABILITY_TYPE`` 
:Default: ``EXPOSURE``
:Type: *keyword*

   - ``EXPOSURE`` {???}
   - ``SEQUENCE`` {???}

-------------------------

:Parameter: ``VERBOSE_TYPE``
:Default: ``NORMAL``
:Type: *keyword*

   Degree of verbosity of the software on screen:

   - ``QUIET`` No Output besides warnings and error messages
   - ``NORMAL`` ``Normal`` display with messages updated in real time using
     ASCII escapes-sequences
   - ``LOG`` Like ``NORMAL``, but without real-time messages and
     ASCII escape-sequences
   - ``FULL`` Everything

-------------------------

:Parameter: ``WRITE_XML``
:Default: ``Y``
:Type: *Boolean*

  If true (``Y``), an XML summary file will be written after completing
  the processing.

-------------------------

:Parameter: ``XML_NAME``
:Default: ``psfex.xml``
:Type: *string*

  File name for the XML output of |PSFEx|.

-------------------------

:Parameter: ``XSL_URL`` 
:Default: ``.``
:Type: *string*

  URL of an XSL style-sheet for the XML output of |PSFEx|. This URL
  will appear in the ``href`` attribute of the ``style-sheet`` tag.


.. include:: keys.rst

