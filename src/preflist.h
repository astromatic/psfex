 /*
 				preflist.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Keywords for the configuration file.
*
*	Last modify:	15/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "key.h"
#ifndef	_CHECK_H_
#include "check.h"
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
  {"CHECKIMAGE_CUBE", P_BOOL, &prefs.check_cubeflag},
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "BASIS", "CHI", "PROTOTYPES", "RESIDUALS", "RAWDATA", "SAMPLES",
	"SNAPSHOTS", "SNAPSHOTS_IMRES", "WEIGHTS",
	"MOFFAT", "-MOFFAT", "-SYMMETRICAL",""},
   0, MAXCHECK, &prefs.ncheck_type},
  {"HOMOBASIS_NUMBER", P_INT, &prefs.homobasis_number, 0,10000},
  {"HOMOBASIS_SCALE", P_FLOAT, &prefs.homobasis_scale, 0,0, 0.0,1.0e3},
  {"HOMOBASIS_TYPE", P_KEY, &prefs.homobasis_type, 0,0, 0.0,0.0,
   {"NONE", "GAUSS-LAGUERRE", ""}},
  {"HOMOKERNEL_NAME", P_STRING, prefs.homokernel_name},
  {"HOMOPSF_PARAMS", P_FLOATLIST, prefs.homopsf_params, 0,0, 0.0,100.0, {""},
     2,2, &prefs.nhomopsf_params},
  {"NEWBASIS_TYPE", P_KEY, &prefs.newbasis_type, 0,0, 0.0,0.0,
	{"NONE", "PCA_MULTI", "PCA_SINGLE", ""}},
  {"NEWBASIS_NUMBER", P_INT, &prefs.newbasis_number, 0,1000},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"PSFVAR_DEGREES", P_INTLIST, prefs.group_deg, 0,32,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ngroup_deg},
  {"PSFVAR_KEYS", P_STRINGLIST, prefs.context_name, 0,0,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_name},
  {"PSFVAR_GROUPS", P_INTLIST, prefs.context_group, 1,MAXCONTEXT,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_group},
  {"PSFVAR_NSNAP", P_INT, &prefs.context_nsnap, 1,16},
  {"PSF_ACCURACY", P_FLOAT, &prefs.prof_accuracy, 0,0, 0.0,1.0},
  {"PSF_NAME", P_STRING, prefs.psf_name},
  {"PSF_RECENTER", P_BOOL, &prefs.recenter_flag},
  {"PSF_SAMPLING", P_FLOAT, &prefs.psf_step, 0,0, 0.0,1.0e3},
  {"PSF_SIZE", P_INTLIST, prefs.retisize, 1,1024, 0.0,0.0, {""},
     1,2, &prefs.nretisize},
  {"SAMPLE_AUTOSELECT", P_BOOL, &prefs.autoselect_flag},
  {"SAMPLE_FLAGMASK", P_INT, &prefs.flag_mask, 0,0xffff},
  {"SAMPLE_FWHMRANGE", P_FLOATLIST, prefs.fwhmrange, 0,0, 0.0,1e3, {""},
     2,2, &prefs.nfwhmrange},
  {"SAMPLE_MAXELONG", P_FLOAT, &prefs.maxelong, 0,0, 1.0, BIG},
  {"SAMPLE_MINSN", P_FLOAT, &prefs.minsn, 0,0, 1e-6,1e15},
  {"SAMPLE_VARIABILITY", P_FLOAT, &prefs.maxvar, 0,0, 0.0, BIG},
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
"PSF_NAME        default.psf     # Output PSF filename",
"BASIS_TYPE      PIXEL_AUTO      # NONE, PIXEL, GAUSS-LAGUERRE or FILE",
"BASIS_NUMBER    16              # Basis number or parameter",
"*BASIS_NAME      basis.fits      # Basis filename (FITS data-cube)",
"*BASIS_SCALE     1.0             # Gauss-Laguerre beta parameter",
"*NEWBASIS_TYPE   NONE            # Create new basis: NONE,PCA_MULTI or PCA_SINGLE",
"*NEWBASIS_NUMBER 8               # Number of new basis vectors",
"PSF_SAMPLING    0.0             # Sampling step in pixel units (0.0 = auto)",
"PSF_ACCURACY    0.01            # Accuracy to expect from PSF \"pixel\" values",
"PSF_SIZE        25,25           # Image size of the PSF model",
"PSF_RECENTER    Y               # Allow recentering of PSF-candidates Y/N ?",
" ",
"#----------------------------- PSF variability -----------------------------",
" ",
"PSFVAR_KEYS    X_IMAGE,Y_IMAGE # SExtractor or FITS (preceded by :) params",
"PSFVAR_GROUPS  1,1             # Group tag for each context key",
"PSFVAR_DEGREES   1             # Polynom degree for each group",
"*PSFVAR_NSNAP   7               # Number of PSF snapshots per axis",
" ",
"#----------------------------- Sample selection ------------------------------",
" ",
"SAMPLE_AUTOSELECT  Y            # Automatically select the FWHM (Y/N) ?",
"SAMPLE_FWHMRANGE   2.0,10.0     # Allowed FWHM range",
"SAMPLE_VARIABILITY 0.2          # Allowed FWHM variability (1.0 = 100%)",
"SAMPLE_MINSN       20           # Minimum S/N for a source to be used",
"SAMPLE_MAXELONG    2.0          # Maximum A/B for a source to be used",
"*SAMPLE_FLAGMASK    0x00fe       # Rejection mask on SExtractor FLAGS",
"*BADPIXEL_FILTER    N            # Filter bad-pixels in samples (Y/N) ?",
"*BADPIXEL_NMAX      0            # Maximum number of bad pixels allowed",
" ",
"*#----------------------- PSF homogeneisation kernel --------------------------",
"*",
"*HOMOBASIS_TYPE     NONE         # NONE or GAUSS-LAGUERRE",
"*HOMOBASIS_NUMBER   10           # Kernel basis number or parameter",
"*HOMOBASIS_SCALE    1.0          # GAUSS-LAGUERRE beta parameter",
"*HOMOPSF_PARAMS     2.0, 3.0     # Moffat parameters of the idealised PSF",
"*HOMOKERNEL_NAME    homo.fits    # Output PSF homogenisation kernel filename",
"*",
"#------------------------------ Check-Images ---------------------------------",
" ",
"CHECKIMAGE_TYPE PROTOTYPES,SAMPLES,RESIDUALS,RAWDATA,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL",
"                                # Check-image types",
"CHECKIMAGE_NAME proto.fits,samp.fits,resi.fits,raw.fits,snap.fits,moffat.fits,submoffat.fits,subsym.fits",
"                                # Check-image filenames",
"*CHECKIMAGE_CUBE N               # Save check-images as datacubes (Y/N) ?",
" ",
"#----------------------------- Miscellaneous ---------------------------------",
" ",
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

