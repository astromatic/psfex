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
*	Last modify:	07/10/2007
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
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "BASIS", "CHI", "PROTOTYPES", "RESIDUALS", "RAWDATA", "SAMPLES",
	"SNAPSHOTS", "SNAPSHOTS_IMRES", "WEIGHTS", "PC_CONVOLVED",
	"MOFFAT", "-MOFFAT", "-SYMMETRICAL",""},
   0, MAXCHECK, &prefs.ncheck_type},
  {"CONTEXT_KEYS", P_STRINGLIST, prefs.context_name, 0,0,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_name},
  {"CONTEXT_GROUPS", P_INTLIST, prefs.context_group, 1,MAXCONTEXT,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_group},
  {"CONTEXT_NSNAP", P_INT, &prefs.context_nsnap, 1,16},
  {"GROUP_DEGREES", P_INTLIST, prefs.group_deg, 0,32,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ngroup_deg},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"PC_INCLUDE", P_BOOL, &prefs.pc_flag},
  {"PC_NAME", P_STRING, prefs.pc_name},
  {"PC_NPC",  P_INT, &prefs.pc_npc, 0,1000000},
  {"PSF_ACCURACY", P_FLOAT, &prefs.prof_accuracy, 0,0, 0.0,1.0},
  {"PSF_AUTOSELECT", P_BOOL, &prefs.autoselect_flag},
  {"PSF_NSUPER", P_INT, &prefs.nsuper,0,1000000000},
  {"PSF_FLAGMASK", P_INT, &prefs.flag_mask, 0,0xffff},
  {"PSF_FWHMRANGE", P_FLOATLIST, prefs.fwhmrange, 0,0, 0.0,1e3, {""},
     2,2, &prefs.nfwhmrange},
  {"PSF_MAXELONG", P_FLOAT, &prefs.maxelong, 0,0, 1.0, BIG},
  {"PSF_MINSN", P_FLOAT, &prefs.minsn, 0,0, 1e-6,1e15},
  {"PSF_NAME", P_STRING, prefs.psf_name},
  {"PSF_RECENTER", P_BOOL, &prefs.recenter_flag},
  {"PSF_SAMPLING", P_FLOAT, &prefs.psf_step, 0,0, 0.0,1.0e3},
  {"PSF_SIZE", P_INTLIST, prefs.retisize, 1,1024, 0.0,0.0, {""},
     1,2, &prefs.nretisize},
  {"PSF_TYPE", P_KEY, &prefs.psf_type, 0,0, 0.0,0.0,
   {"SINC","GAUSS-LAGUERRE",""}},
{"PSF_BETA", P_FLOAT, &prefs.psf_beta, 0,0, 0.0,1.0e3},
  {"PSF_VARIABILITY", P_FLOAT, &prefs.maxvar, 0,0, 0.0, BIG},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","LOG","FULL",""}},
  {"XML_NAME", P_STRING, prefs.xml_name},
  {"XSL_URL", P_STRING, prefs.xsl_name},
  {"WRITE_XML", P_BOOL, &prefs.xml_flag},
  {""}
 };

char		keylist[sizeof(key)/sizeof(pkeystruct)][16];
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
"PSF_TYPE        SINC            # SINC or GAUSS-LAGUERRE",
"*PSF_BETA        1.0             # Gauss-Laguerre beta parameter",
"PSF_NSUPER      16              # Super-resolution number parameter",
"PSF_SAMPLING    0.0             # Sampling step in pixel units (0.0 = auto)",
"PSF_ACCURACY    0.01            # Accuracy to expect from PSF \"pixel\" values",
"PSF_SIZE        25,25           # Image size of the PSF model",
"PSF_RECENTER    Y               # Allow recentering of PSF-candidates Y/N ?",
" ",
"#----------------------------- Sample selection ------------------------------",
" ",
"PSF_AUTOSELECT  Y               # Automatically select the FWHM (Y/N) ?",
"PSF_FWHMRANGE   2.0,10.0        # Allowed FWHM range",
"PSF_VARIABILITY 0.2             # Allowed PSF variability (1.0 = 100%)",
"PSF_MINSN       20              # Minimum S/N for a source to be used",
"PSF_MAXELONG    2.0             # Maximum A/B for a source to be used",
"*PSF_FLAGMASK    0x00fe          # Rejection mask on SExtractor FLAGS",
"*BADPIXEL_FILTER N               # Filter bad-pixels in samples (Y/N) ?",
"*BADPIXEL_NMAX   0               # Maximum number of bad pixels allowed",
" ",
"#----------------------------- PSF dependencies ------------------------------",
" ",
"CONTEXT_KEYS    X_IMAGE,Y_IMAGE # SExtractor or FITS (preceded by :) params",
"CONTEXT_GROUPS  1,1             # Group tag for each context key",
"GROUP_DEGREES   1               # Polynom degree for each group",
"*CONTEXT_NSNAP   7               # Number of PSF snapshots per axis",
" ",
"#------------------------------ Check-Images ---------------------------------",
" ",
"CHECKIMAGE_TYPE PROTOTYPES,SAMPLES,RESIDUALS,RAWDATA,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL",
"                                # Check-image types",
"CHECKIMAGE_NAME proto.fits,samp.fits,resi.fits,raw.fits,snap.fits,moffat.fits,submoffat.fits,subsym.fits",
"                                # Check-image filenames",
"* ",
"*#---------------------- Galaxy Principal Components --------------------------",
"* ",
"*PC_INCLUDE      N               # Process galaxy Principal Components (Y/N) ?",
"*PC_NAME         default.pc      # File to store galaxy Principal Components",
"*PC_NPC          0               # Number of principal components",
" ",
"#----------------------------- Miscellaneous ---------------------------------",
" ",
"VERBOSE_TYPE    NORMAL          # can be QUIET,NORMAL,LOG or FULL",
"WRITE_XML       Y               # Write XML file (Y/N)?",
"XML_NAME        psfex.xml       # Filename for XML output",
"*XSL_URL        " XSL_URL,
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

