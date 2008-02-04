 /*
 				makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Main program.
*
*	Last modify:	04/02/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"check.h"
#include	"diagnostic.h"
#include	"homo.h"
#include	"pca.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"xml.h"

psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, int diagflag);
time_t		thetime, thetime2;

/********************************** makeit ***********************************/
/*
*/
void	makeit(char **incatnames, int ncat)

  {
   psfstruct		**psf, **cpsf;
   setstruct		*set;
   catstruct		*cat;
   tabstruct		*tab;
   struct tm		*tm;
   static char		str[MAXCHAR];
   float		*psfsteps, *basis,
			psfstep;
   int			c,i, ntab, ext, next, nmed, nbasis;

/* Install error logging */
  error_installfunc(write_error);

/* Processing start date and time */
  thetime = time(NULL);
  tm = localtime(&thetime);
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT,
	"----- %s %s started on %s at %s with %d thread%s\n\n",
		BANNER,
		MYVERSION,
		prefs.sdate_start,
		prefs.stime_start,
		prefs.nthreads,
		prefs.nthreads>1? "s":"");

/* Compute the number of valid input extensions */
  if (!(cat = read_cat(incatnames[0])))
    error(EXIT_FAILURE, "*Error*: cannot open ", incatnames[0]);
  tab = cat->tab;
  next = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
    {
/*--  Check for the next valid image extension */
    if ((tab->naxis != 2)
        || (strncmp(tab->xtension, "BINTABLE", 8)
		&& strncmp(tab->xtension, "ASCTABLE", 8))
	|| (strncmp(tab->extname, "LDAC_OBJECTS", 8)
		&& strncmp(tab->extname, "OBJECTS", 8)))
      continue;
    next++;
    }
  free_cat(&cat, 1);
  if (prefs.xml_flag)
    init_xml(next);  

  psfstep = prefs.psf_step;
  psfsteps = NULL;
  nbasis = 0;
  basis = NULL;

  if (prefs.newbasis_type==NEWBASIS_PCAMULTI && (!psfstep))
    {
/* A first run through all samples only to derive a common pixel step */
    QMALLOC(psfsteps, float, next);
    for (ext=0 ; ext<next; ext++)
      {
      set = load_samples(incatnames, ncat, ext, next);
      psfsteps[ext] = (float)(psfstep? psfstep : (set->fwhm/2.35)*0.5);
      end_set(set);
      }
    }
  else if (prefs.newbasis_type==NEWBASIS_PCASINGLE)
    {
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT,
	"----- First pass for Principal Component Analysis (single): "
	"%dx%d PSFs required\n", ncat, next);
    if (!prefs.psf_step)
      {
      QMALLOC(psfsteps, float, next);
      for (ext=0 ; ext<next; ext++)
        {
        set = load_samples(incatnames, ncat, ext, next);
        psfsteps[ext] = (float)((set->fwhm/2.35)*0.5);
        end_set(set);
        }
      psfstep = fast_median(psfsteps, next);
      }
/*-- Derive a new common PCA basis for all extensions */
    QMALLOC(cpsf, psfstruct *, ncat*next);
    for (ext=0 ; ext<next; ext++)
      for (c=0; c<ncat; c++)
        {
        set = load_samples(&incatnames[c], 1, ext, next);
        sprintf(str, "Computing PSF model for catalog %s...", incatnames[c]);
        NFPRINTF(OUTPUT, str);
        cpsf[c+ext*ncat] = make_psf(set, psfstep, NULL, 0, 0);
        end_set(set);
        }
    nbasis = prefs.newbasis_number;
    basis = pca_make(cpsf, ncat*next, nbasis);
    for (i=0 ; i<ncat*next; i++)
      psf_end(cpsf[i]);
    free(cpsf);
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT,
	"----- Second pass for Principal Component Analysis (single):\n\n");
    }

/* Create an array of PSFs (one PSF for each extension) */
  QMALLOC(psf, psfstruct *, next);
  QIPRINTF(OUTPUT,
        " extension accepted/total sampling chi2/dof FWHM(pix) elong."
	" residuals asymmetry");
  for (ext=0 ; ext<next; ext++)
    {
    if (prefs.newbasis_type == NEWBASIS_PCAMULTI)
/*---- Derive a new PCA basis for each extension */
      {
      psfstep = prefs.psf_step? prefs.psf_step : psfsteps[ext];
      QMALLOC(cpsf, psfstruct *, ncat);
      for (c=0; c<ncat; c++)
        {
        set = load_samples(&incatnames[c], 1, ext, next);
        sprintf(str, "Computing PSF model for catalog %s...", incatnames[c]);
        NFPRINTF(OUTPUT, str);
        cpsf[c] = make_psf(set, psfstep, NULL, 0, 0);
        end_set(set);
        }
      nbasis = prefs.newbasis_number;
      basis = pca_make(cpsf, ncat, nbasis);
      for (c=0 ; c<ncat; c++)
        psf_end(cpsf[c]);
      free(cpsf);
      }

/*-- Load all the samples */
    set = load_samples(incatnames, ncat, ext, next);
    if (prefs.newbasis_type == NEWBASIS_NONE && !psfstep)
      psfstep = (float)((set->fwhm/2.35)*0.5);
    psf[ext] = make_psf(set, psfstep, basis, nbasis, 1);
    nmed = psf[ext]->nmed;
    QPRINTF(OUTPUT, "[%3d/%-3d]     %5d/%-5d  %6.2f    %6.2f %6.2f    %5.3f"
	"    %5.2f     %5.2f\n",
	ext+1, next,
	psf[ext]->samples_accepted, psf[ext]->samples_loaded,
	psf[ext]->pixstep,
	psf[ext]->chi2,
	sqrt(psf[ext]->moffat[nmed].fwhm_min*psf[ext]->moffat[nmed].fwhm_max),
	psf[ext]->moffat[nmed].fwhm_max/psf[ext]->moffat[nmed].fwhm_min,
	psf[ext]->moffat[nmed].residuals, psf[ext]->moffat[nmed].symresiduals);
/*-- Save result */
    NFPRINTF(OUTPUT,"Saving the PSF description...");
    psf_save(psf[ext], prefs.psf_name, ext, next);

/* Create a homogenisation kernel */
    if (prefs.homobasis_type != HOMOBASIS_NONE)
      {
      NFPRINTF(OUTPUT, "Computing homogenisation kernel...");
      psf_homo(psf[ext], prefs.homokernel_name, prefs.homopsf_params,
		prefs.homobasis_number, prefs.homobasis_scale, ext, next);
      }

/*-- Save "Check-images" */
    for (i=0; i<prefs.ncheck_type; i++)
      if (prefs.check_type[i])
        {
        sprintf(str, "Saving CHECK-image #%d...", i+1);
        NFPRINTF(OUTPUT, str);
        psf_writecheck(psf[ext], set, prefs.check_name[i],
		prefs.check_type[i], ext, next, prefs.check_cubeflag);
        }
/*-- Update XML */
    if (prefs.xml_flag)
      update_xml(psf[ext], ncat);
/*-- Free memory */
    end_set(set);
    }

  free(psfsteps);

/* Processing end date and time */
  thetime2 = time(NULL);
  tm = localtime(&thetime2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = difftime(thetime2, thetime);

/* Write XML */
  if (prefs.xml_flag)
    {
    NFPRINTF(OUTPUT, "Writing XML file...");
    write_xml(prefs.xml_name);
    end_xml();
    }

/* Free memory */
  for (ext = 0 ; ext<next; ext++)
    psf_end(psf[ext]);
  free(psf);

  return;
  }


/****** make_psf *************************************************************
PROTO	psfstruct	*make_psf(setstruct *set, float psfstep,
				float *basis, int nbasis, int diagflag)
PURPOSE	Make PSFs from a set of FITS binary catalogs.
INPUT	Pointer to a sample set,
	PSF sampling step,
	Pointer to basis image vectors,
	Number of basis vectors,
	Diagnostic flag.
OUTPUT  Pointer to the PSF structure.
NOTES   Diagnostics are computed only if diagflag != 0.
AUTHOR  E. Bertin (IAP)
VERSION 15/01/2008
 ***/
psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, int diagflag)
  {
   psfstruct		*psf;
   basistypenum		basistype;

  NFPRINTF(OUTPUT,"Initializing PSF modules...");
  psf = psf_init(prefs.context_name, prefs.context_group,
		prefs.ncontext_group,
		prefs.group_deg, prefs.ngroup_deg,
		set->retisize[0], set->retisize[1],
		psfstep,
		set->nsample);

  psf->samples_loaded = set->nsample;
  psf->fwhm = set->fwhm;
  
/* Make the basic PSF-model (1st pass) */
  NFPRINTF(OUTPUT,"Modeling the PSF.");
  psf_make(psf, set);
  if (basis && nbasis)
    {
    QMEMCPY(basis, psf->basis, float, nbasis*psf->size[0]*psf->size[1]);
    psf->nbasis = nbasis;
    }
  else
    {
    NFPRINTF(OUTPUT,"Generating the PSF model...");
    basistype = prefs.basis_type;
    if (basistype==BASIS_PIXEL_AUTO)
      basistype = (psf->fwhm < PSF_AUTO_FWHM)? BASIS_PIXEL : BASIS_NONE;
    psf_makebasis(psf, set, basistype, prefs.basis_number);
    }
  psf_refine(psf, set);

/* Remove bad PSF candidates */
  if (set->nsample>1)
    {
    psf_clean(psf, set);

/*-- Make the basic PSF-model (2nd pass) */
    NFPRINTF(OUTPUT,"Modeling the PSF.");
    psf_make(psf, set);
    psf_refine(psf, set);
    }

/* Remove bad PSF candidates */
  if (set->nsample>1)
    psf_clean(psf, set);

  psf->samples_accepted = set->nsample;

/* Refine the PSF-model */
  psf_refine(psf, set);

/* Clip the PSF-model */
  psf_clip(psf);

/*-- Just check the Chi2 */
  psf->chi2 = set->nsample? psf_chi2(psf, set) : 0.0;

/* Make a diagnostic of the PSF */
  if (diagflag)
    psf_diagnostic(psf);

  return psf;
  }


/****** write_error ********************************************************
PROTO	void    write_error(char *msg1, char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/02/2007
 ***/
void	write_error(char *msg1, char *msg2)
  {
   char	error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    write_xmlerror(prefs.xml_name, error);
  end_xml();

  return;
  }


