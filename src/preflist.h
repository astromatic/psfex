
/*
*				preflist.h
*
* Configuration keyword definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		26/02/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "key.h"
#ifndef	_CHECK_H_
#include "check.h"
#endif

#ifndef	_CONTEXT_H_
#include "context.h"
#endif

#ifndef _XML_H_
#include "xml.h"
#endif

#ifdef  USE_THREADS
#define	THREADS_PREFMAX	THREADS_NMAX
#else
#define	THREADS_PREFMAX	65535
#endif

/*-------------------------------- initialization ---------------------------*/
int	idummy;

pkeystruct key[] =
 {
  {"BADPIXEL_FILTER", P_BOOL, &prefs.badpix_flag},
  {"BADPIXEL_NMAX", P_INT, &prefs.badpix_nmax, 0,100000000},
  {"BASIS_NAME", P_STRING, prefs.basis_name},
  {"BASIS_NUMBER", P_INT, &prefs.basis_number, 0,10000},
  {"BASIS_SCALE", P_FLOAT, &prefs.basis_scale, 0,0, 0.0,1.0e3},
  {"BASIS_TYPE", P_KEY, &prefs.basis_type, 0,0, 0.0,0.0,
   {"NONE", "PIXEL", "GAUSS-LAGUERRE", "FILE", "PIXEL_AUTO", ""}},
  {"CENTER_KEYS", P_STRINGLIST, prefs.center_key, 0,0,0.0,0.0,
    {""}, 2, 2, &prefs.ncenter_key},
  {"CHECKIMAGE_CUBE", P_BOOL, &prefs.check_cubeflag},
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "BASIS", "CHI", "PROTOTYPES", "RESIDUALS", "RESIDUALS_GRID",
	"SAMPLES", "SAMPLES_GRID",
	"SNAPSHOTS", "SNAPSHOTS_IMRES", "WEIGHTS",
	"MOFFAT", "-MOFFAT", "-SYMMETRICAL", "GREAT10", ""},
   0, MAXCHECK, &prefs.ncheck_type},
  {"CHECKPLOT_ANTIALIAS", P_BOOL, &prefs.cplot_antialiasflag},
  {"CHECKPLOT_DEV", P_KEYLIST, prefs.cplot_device, 0,0, 0.0,0.0,
    {"NULL", "XWIN", "TK", "XTERM", "PLMETA", "PS", "PSC", "XFIG", "LJIIP",
	"LJ_HPGL", "IMP", "PBM", "PNG", "JPEG", "PSTEX", "AQT", "PDF", "SVG",
	""},
    0, MAXCHECK, &prefs.ncplot_device},
  {"CHECKPLOT_NAME", P_STRINGLIST, prefs.cplot_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncplot_name},
  {"CHECKPLOT_RES", P_INTLIST, prefs.cplot_res, 0,16384,0.0,0.0,
    {""}, 1, 2, &prefs.ncplot_res},
   {"CHECKPLOT_TYPE", P_KEYLIST, prefs.cplot_type, 0,0, 0.0,0.0,
    {"NONE", "FWHM", "ELLIPTICITY", "MOFFAT_RESIDUALS", "ASYMMETRY",
	"COUNTS", "COUNT_FRACTION", "CHI2", "RESIDUALS", "GREAT10", ""},
    0, MAXCHECK, &prefs.ncplot_type},
  {"HIDDENMEF_TYPE", P_KEY, &prefs.hidden_mef_type, 0,0, 0.0,0.0,
	{"INDEPENDENT", "COMMON", ""}},
  {"HOMOBASIS_NUMBER", P_INT, &prefs.homobasis_number, 0,10000},
  {"HOMOBASIS_SCALE", P_FLOAT, &prefs.homobasis_scale, 0,0, 0.0,1.0e3},
  {"HOMOBASIS_TYPE", P_KEY, &prefs.homobasis_type, 0,0, 0.0,0.0,
   {"NONE", "GAUSS-LAGUERRE", ""}},
  {"HOMOKERNEL_DIR", P_STRING, prefs.homokernel_dir},
  {"HOMOKERNEL_SUFFIX", P_STRING, prefs.homokernel_suffix},
  {"HOMOPSF_PARAMS", P_FLOATLIST, prefs.homopsf_params, 0,0, 0.0,100.0, {""},
     2,2, &prefs.nhomopsf_params},
  {"MEF_TYPE", P_KEY, &prefs.psf_mef_type, 0,0, 0.0,0.0,
	{"INDEPENDENT", "COMMON", ""}},
  {"NEWBASIS_TYPE", P_KEY, &prefs.newbasis_type, 0,0, 0.0,0.0,
	{"NONE", "PCA_INDEPENDENT", "PCA_COMMON", ""}},
  {"NEWBASIS_NUMBER", P_INT, &prefs.newbasis_number, 0,1000},
  {"NTHREADS", P_INT, &prefs.nthreads, -THREADS_PREFMAX, THREADS_PREFMAX},
  {"OUTCAT_NAME", P_STRING, prefs.outcat_name},
  {"OUTCAT_TYPE", P_KEY, &prefs.outcat_type, 0,0, 0.0,0.0,
   {"NONE", "ASCII_HEAD", "ASCII", "ASCII_VOTABLE", "FITS_LDAC", ""}},
  {"PHOTFLUX_KEY", P_STRING, prefs.photflux_key},
  {"PHOTFLUXERR_KEY", P_STRING, prefs.photfluxerr_key},
  {"PSFVAR_DEGREES", P_INTLIST, prefs.group_deg, 0,32,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ngroup_deg},
  {"PSFVAR_KEYS", P_STRINGLIST, prefs.context_name, 0,0,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_name},
  {"PSFVAR_GROUPS", P_INTLIST, prefs.context_group, 1,MAXCONTEXT,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_group},
  {"PSFVAR_NSNAP", P_INT, &prefs.context_nsnap, 1,256},
  {"PSF_ACCURACY", P_FLOAT, &prefs.prof_accuracy, 0,0, 0.0,1.0},
  {"PSF_DIR", P_STRING, prefs.psf_dir},
  {"PSF_PIXELSIZE", P_FLOATLIST, prefs.psf_pixsize, 0,0, 0.0,100.0, {""},
     1,2, &prefs.npsf_pixsize},
  {"PSF_RECENTER", P_BOOL, &prefs.recenter_flag},
  {"PSF_SAMPLING", P_FLOAT, &prefs.psf_step, 0,0, 0.0,1.0e3},
  {"PSF_SIZE", P_INTLIST, prefs.psf_size, 1,1024, 0.0,0.0, {""},
     1,2, &prefs.npsf_size},
  {"PSF_SUFFIX", P_STRING, prefs.psf_suffix},
  {"SAMPLE_AUTOSELECT", P_BOOL, &prefs.autoselect_flag},
  {"SAMPLE_FLAGMASK", P_INT, &prefs.flag_mask, 0,0xffff},
  {"SAMPLE_FWHMRANGE", P_FLOATLIST, prefs.fwhmrange, 0,0, 0.0,1e3, {""},
     2,2, &prefs.nfwhmrange},
  {"SAMPLE_IMAFLAGMASK", P_INT, &prefs.imaflag_mask, 0,0xff, 0.0,0.0},
  {"SAMPLE_MAXELLIP", P_FLOAT, &prefs.maxellip, 0,0, 0.0, 1.0},
  {"SAMPLE_MINSN", P_FLOAT, &prefs.minsn, 0,0, 1e-6,1e15},
//  {"SAMPLE_NMAX", P_INT, &prefs.nmax, 0,2147483648},
  {"SAMPLE_VARIABILITY", P_FLOAT, &prefs.maxvar, 0,0, 0.0, BIG},
  {"SAMPLE_WFLAGMASK", P_INT, &prefs.wflag_mask, 0,0xff, 0.0,0.0},
  {"SAMPLEVAR_TYPE", P_KEY, &prefs.var_type, 0,0, 0.0,0.0,
	{"NONE", "SEEING",""}},
  {"STABILITY_TYPE", P_KEY, &prefs.stability_type, 0,0, 0.0,0.0,
	{"EXPOSURE", "SEQUENCE", ""}},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","LOG","FULL",""}},
  {"XML_NAME", P_STRING, prefs.xml_name},
  {"XSL_URL", P_STRING, prefs.xsl_name},
  {"WRITE_XML", P_BOOL, &prefs.xml_flag},
  {""}
 };

char		keylist[sizeof(key)/sizeof(pkeystruct)][32];
const char	notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
" ",
"#-------------------------------- PSF model ----------------------------------",
" ",
"BASIS_TYPE      PIXEL_AUTO      # NONE, PIXEL, GAUSS-LAGUERRE or FILE",
"BASIS_NUMBER    20              # Basis number or parameter",
"*BASIS_NAME      basis.fits      # Basis filename (FITS data-cube)",
"*BASIS_SCALE     1.0             # Gauss-Laguerre beta parameter",
"*NEWBASIS_TYPE   NONE            # Create new basis: NONE, PCA_INDEPENDENT",
"*                                # or PCA_COMMON",
"*NEWBASIS_NUMBER 8               # Number of new basis vectors",
"PSF_SAMPLING    0.0             # Sampling step in pixel units (0.0 = auto)",
"*PSF_PIXELSIZE   1.0             # Effective pixel size in pixel step units",
"PSF_ACCURACY    0.01            # Accuracy to expect from PSF \"pixel\" values",
"PSF_SIZE        25,25           # Image size of the PSF model",
"*PSF_RECENTER    N               # Allow recentering of PSF-candidates Y/N ?",
"*MEF_TYPE        INDEPENDENT     # INDEPENDENT or COMMON",
" ",
"#------------------------- Point source measurements -------------------------",
" ",
"CENTER_KEYS     X_IMAGE,Y_IMAGE # Catalogue parameters for source pre-centering",
"PHOTFLUX_KEY    FLUX_APER(1)    # Catalogue parameter for photometric norm.",
"PHOTFLUXERR_KEY FLUXERR_APER(1) # Catalogue parameter for photometric error",
" ",
"#----------------------------- PSF variability -------------------------------",
" ",
"PSFVAR_KEYS     X_IMAGE,Y_IMAGE # Catalogue or FITS (preceded by :) params",
"PSFVAR_GROUPS   1,1             # Group tag for each context key",
"PSFVAR_DEGREES  2               # Polynom degree for each group",
"*PSFVAR_NSNAP    9               # Number of PSF snapshots per axis",
"*HIDDENMEF_TYPE  COMMON          # INDEPENDENT or COMMON",
"*STABILITY_TYPE  EXPOSURE        # EXPOSURE or SEQUENCE",
" ",
"#----------------------------- Sample selection ------------------------------",
" ",
"SAMPLE_AUTOSELECT  Y            # Automatically select the FWHM (Y/N) ?",
"SAMPLEVAR_TYPE     SEEING       # File-to-file PSF variability: NONE or SEEING",
"SAMPLE_FWHMRANGE   2.0,10.0     # Allowed FWHM range",
"SAMPLE_VARIABILITY 0.2          # Allowed FWHM variability (1.0 = 100%)",
"SAMPLE_MINSN       20           # Minimum S/N for a source to be used",
"SAMPLE_MAXELLIP    0.3          # Maximum (A-B)/(A+B) for a source to be used",
"*SAMPLE_FLAGMASK    0x00fe       # Rejection mask on SExtractor FLAGS",
"*SAMPLE_WFLAGMASK   0x0000       # Rejection mask on SExtractor FLAGS_WEIGHT",
"*SAMPLE_IMAFLAGMASK 0x0          # Rejection mask on SExtractor IMAFLAGS_ISO",
//"*SAMPLE_NMAX        0            # Maximum number of samples per extension",
"*BADPIXEL_FILTER    N            # Filter bad-pixels in samples (Y/N) ?",
"*BADPIXEL_NMAX      0            # Maximum number of bad pixels allowed",
" ",
"*#----------------------- PSF homogeneisation kernel --------------------------",
"*",
"*HOMOBASIS_TYPE     NONE         # NONE or GAUSS-LAGUERRE",
"*HOMOBASIS_NUMBER   10           # Kernel basis number or parameter",
"*HOMOBASIS_SCALE    1.0          # GAUSS-LAGUERRE beta parameter",
"*HOMOPSF_PARAMS     2.0, 3.0     # Moffat parameters of the idealised PSF",
"*HOMOKERNEL_DIR                  # Where to write kernels (empty=same as input)",
"*HOMOKERNEL_SUFFIX  .homo.fits   # Filename extension for homogenisation kernels",
"*",
"*#----------------------------- Output catalogs -------------------------------",
"*",
"*OUTCAT_TYPE        NONE         # NONE, ASCII_HEAD, ASCII, FITS_LDAC",
"*OUTCAT_NAME        psfex_out.cat  # Output catalog filename",
"*",
"#------------------------------- Check-plots ----------------------------------",
" ",
"CHECKPLOT_DEV       PNG         # NULL, XWIN, TK, PS, PSC, XFIG, PNG,",
"                                # JPEG, AQT, PDF or SVG",
"*CHECKPLOT_RES       0           # Check-plot resolution (0 = default)",
"*CHECKPLOT_ANTIALIAS Y           # Anti-aliasing using convert (Y/N) ?",
"CHECKPLOT_TYPE      FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS",
"                                # or NONE",
"CHECKPLOT_NAME      fwhm, ellipticity, counts, countfrac, chi2, resi",
" ",
"#------------------------------ Check-Images ---------------------------------",
" ",
"CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS",
"                                # or MOFFAT,-MOFFAT,-SYMMETRICAL",
"CHECKIMAGE_NAME chi.fits,proto.fits,samp.fits,resi.fits,snap.fits",
"                                # Check-image filenames",
"*CHECKIMAGE_CUBE N               # Save check-images as datacubes (Y/N) ?",
" ",
"#----------------------------- Miscellaneous ---------------------------------",
" ",
"PSF_DIR                         # Where to write PSFs (empty=same as input)",
"*PSF_SUFFIX      .psf            # Filename extension for output PSF filename",
"VERBOSE_TYPE    NORMAL          # can be QUIET,NORMAL,LOG or FULL",
"WRITE_XML       Y               # Write XML file (Y/N)?",
"XML_NAME        psfex.xml       # Filename for XML output",
"*XSL_URL         " XSL_URL,
"*                                # Filename for XSL style-sheet",
#ifdef USE_THREADS
"NTHREADS        0               # Number of simultaneous threads for",
"                                # the SMP version of " BANNER,
"                                # 0 = automatic",
#else
"NTHREADS        1               # 1 single thread",
#endif
" ",
""
 };

