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


License
=======

|PSFEx|_ *is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.* |PSFEx| *is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details. You should have received a copy
of the GNU General Public License along with* |PSFEx|. *If not, see*
`www.gnu.org/licenses/ <http://www.gnu.org/licences/>`_.

Installing the software
=======================

Hardware requirements
---------------------

|PSFEx| runs in (ANSI) text-mode from a shell. A window system is
not necessary for basic operation.

When it comes to memory usage, the amount required by |PSFEx| depends mostly on
the number of point sources present in the input catalogues times the number of
pixels in the small image that represents each of them. A typical figure is
about 15 kbytes per point source; hence even on a modest computer with
1GB of memory, more than 20,000 point sources can easily be
accommodated at once.

Note that |PSFEx| takes advantage of multiple CPU cores for some operations.

Obtaining |PSFEx|
-----------------

For Linux users, the simplest way to have |PSFEx| up and running is to install
the standard binary package the comes with your Linux distribution. Run, e.g.,
``apt-get psfex`` (on Debian) or ``dnf psfex`` (Fedora) and
|PSFEx|, as well as all its dependencies, will automatically be installed. If you
decided to install the package this way you may skip the following and move
straight to the :ref:`next section <gettingStarted>`.

However if |PSFEx| is not available in your distribution, or to obtain the most
recent version, the |PSFEx| source package can be downloaded from `the official
GitHub repository
<https://github.com/astromatic/psfex>`_ . One may choose `one of the stable
releases <https://github.com/astromatic/psfex/releases>`_, or for the fearless,
`a copy of the current master development branch
<https://github.com/astromatic/psfex/archive/master.zip>`_.

Software requirements
---------------------

|PSFEx| has been developed on `GNU/Linux <http://en.wikipedia.org/wiki/Linux>`_
machines and should compile on any
`POSIX <http://en.wikipedia.org/wiki/POSIX>`_-compliant system (this includes
|OSX|_ and `Cygwin <http://www.cygwin.com>`_ on |Windows|_, at the price of
some difficulties with the configuration), provided
that the *development* packages of the following libraries have been installed:

* |ATLAS|_ V3.6 and above [#atlas_install]_,
* |FFTw|_ V3.0 and above [#fftw_install]_, 
* |PLPlot|_ V5.9 and above.

On Fedora/Redhat distributions for instance, the development packages above are
available as ``atlas-devel``, ``fftw-devel`` and ``plplot-devel``.
|PLPlot| is only required for producing diagnostic plots. Note
that |ATLAS| and |FFTw| are not necessary if |PSFEx| is linked with
|Intel|'s |MKL|_ library.

Installation
------------

To install from the |GitHub| source package, you must first uncompress the
archive:

.. code-block:: console

  $ unzip psfex-<version>.zip

A new directory called ``psfex-<version>`` should now appear at the current
location on your disk. Enter the directory and generate the files required by
the `autotools <http://en.wikipedia.org/wiki/GNU_Build_System>`_, which the
package relies on:

.. code-block:: console

  $ cd psfex-<version>
  $ sh autogen.sh

A ``configure`` script is created. This script has many options, which may be
listed with the :option:`--help` option:

.. code-block:: console

  $ ./configure --help

No options are required for compiling with the default GNU C compiler (``gcc``)
if all the required libraries are installed at their default locations:

.. code-block:: console

  $ ./configure

Compared to ``gcc`` and the librairies above, the combination of the |Intel|
compiler (``icc``) and the |MKL|_ libraries can give the |PSFEx| executable a
strong boost in performance, thanks to better vectorized code.
If ``icc`` and the |MKL| are installed on your system [#geticc]_ , you can take
advantage of them using

.. code-block:: console

  $ ./configure --enable-mkl

Additionally, if the |PSFEx| binary is to be run on a different machine
that does not have ``icc`` and the |MKL| installed (e.g., a cluster computing
node), you must configure a partially statically linked executable using

.. code-block:: console

  $ ./configure --enable-mkl --enable-auto-flags --enable-best-link

In all cases, |PSFEx| can now be compiled with

.. code-block:: console

  $ make -j

An ``src/psfex`` executable is created. For system-wide installation, run the
usual

.. code-block:: console

  $ sudo make install

You may now check that the software is properly installed by simply
typing in your shell:

.. code-block:: console

  $ psfex

which will return the version number and other basic information (note that some shells require the ``rehash`` command to be run before making a
freshly installed executable accessible in the execution path).

.. [#mac_install] Mac OS X |.dmg|_ packages should be available soon.
.. [#atlas_install] Use the :option:`--with-atlas` and/or
   :option:`--with-atlas-incdir` options of the |PSFEx| :command:`configure`
   script to specify the |ATLAS| library and include paths if |ATLAS| files are 
   installed at unusual locations.
.. [#fftw_install] Make sure that |FFTW| has been compiled with
   :command:`configure` options :option:`--enable-threads --enable-float`.
.. [#geticc] The Linux versions of the |Intel| compiler and |MKL| are
   `available for free to academic researchers, students, educators and open
   source contributors <http://software.intel.com/qualify-for-free-software>`_.

.. include:: keys.rst

