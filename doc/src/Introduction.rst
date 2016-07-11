.. File Introduction.rst

Introduction
===================


What is |PSFEx|?
================

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

|PSFEx| is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version. |PSFEx| is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details. You should have received a copy
of the GNU General Public License along with |PSFEx|. If not, see
`www.gnu.org/licenses/ <http://www.gnu.org/licences/>`_.

Installing the software
=======================

Obtaining |PSFEx|
-----------------

The easiest way to obtain |PSFEx| is to download it from `the official Web site 
<http://astromatic.net/software/psfex>`_ . At this address, the latest versions
of the program (source code, configuration files, and documentation) are
available as standard ``.tar.gz`` Unix source archives as well as |RPM|_ binary
packages for various architectures[#dmg]_. 

Software and hardware requirements
----------------------------------

|PSFEx| has been developed on Unix machines (GNU/Linux), and should compile on
any POSIX-compliant system (this should include Mac OS X and Cygwin under
Windows, at the price of some difficulties with the configuration), provided
that the following libraries/packages have been installed:

* |ATLAS|_ V3.6 and above [#atlas_install]_,
* |FFTw|_ V3.0 and above [#fftw_install]_, 
* |PLPlot|_ V5.9 and above.

|PLPlot| is only required for producing diagnostic plots. Note
that |ATLAS| and |FFTw| are not necessary for the binary
versions of |PSFEx| which come with these libraries statically linked.

The software is run in (ANSI) text-mode from a shell. A window system is
necessary only when |PLPlot| is used in interactive mode.

The amount of memory required by |PSFEx| depends mostly on the number
of point sources present in the input catalogues times the number of pixels
in the small image that represents each of them. A typical figure is
about 15 kbytes per point source; hence even on a modest computer with
256MB of memory, more than 10,000 point sources can easily be
accommodated at once.

Installation from the source archive
------------------------------------

To install from the source, you must first uncompress and “untar” the
archive:

.. code-block:: console

  $ tar zxvf psfex-<version>.tar.gz

A new directory called ``psfex-<version>`` should now appear at the current
location on your disk. You should then enter the directory and follow
the instructions of the ``INSTALL`` file.

Installation from an RPM archive
--------------------------------

|PSFEx| is also available as a binary RPM package for both Linux INTEL x86
(32-bit) and x86-64 (64-bit) architectures. To check which matches your
system, use the shell command:

.. code-block:: console

  $ uname -a

The RPM version of |PSFEx| requires the |PLPlot| package. Make sure it is
installed before proceeding. To install |PSFEx|, type as a root user the
following command in your shell (preceded with su if you don’t have root
access but the system administrator trusts you well enough to make you
part of the wheel group):

.. code-block:: console

  $ rpm -U psfex-<version>-1.<arch>.rpm

It is sometimes necessary to force installation with:

.. code-block:: console

  $ rpm -U --force --nodeps psfex-<version>-1.<arch>.rpm

You may now check that the software is properly installed by simply
typing in your shell:

.. code-block:: console

  $ psfex

(note that some shells require the ``rehash`` command to be run before making a
freshly installed executable accessible in the execution path).

.. [#mac_install] Mac OS X |.dmg|_ packages should be available soon.
.. [#atlas_install] Use the :option:`--with-atlas` and/or
   :option:`--with-atlas-incdir` options of the |PSFEx| :command:`configure`
   script to specify the |ATLAS| library and include paths if |ATLAS| files are 
   installed at unusual locations.
.. [#fftw_install] Make sure that |FFTW| has been compiled with
   :command:`configure` options :option:`--enable-threads --enable-float`.

.. include:: keys.rst

