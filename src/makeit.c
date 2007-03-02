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
*	Last modify:	02/03/2007
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
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"
#include	"xml.h"

time_t		thetime, thetime2;

/********************************** makeit ***********************************/
/*
*/
void	makeit(char **incatnames, int ncat)

  {
   setstruct		*set;
   psfstruct		**psf;
   pcstruct		*pc, *pcc, *pco;
   catstruct		*cat;
   tabstruct		*tab;
   static char		str[MAXCHAR];
   struct tm		*tm;
   float		psfstep;
   int			i, ntab, ext, next, nmed;

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

  QIPRINTF(OUTPUT,
        " extension accepted/total sampling chi2/dof FWHM(pix) elongation"
	" residuals");

/* Create an array of PSFs (one PSF for each extension) */
  QMALLOC(psf, psfstruct *, next);

  for (ext = 0 ; ext<next; ext++)
    {
/*-- Load all the samples */
    NFPRINTF(OUTPUT,"Loading samples...");
    set = load_samples(incatnames, ncat, ext, next);

    sprintf(str, "%d samples loaded.", set->nsample);
    NFPRINTF(OUTPUT, str);

    if (!set->nsample)
      warning("No appropriate source found!!","");

    psfstep = (float)(prefs.psf_step? prefs.psf_step
				: (set->fwhm/2.35)*(1.0-1.0/INTERPFAC));

/*
    NFPRINTF(OUTPUT, "");
    NPRINTF(OUTPUT, "PSF sampled each %.2f pixel%s (%s)\n",
			psfstep,
			psfstep>=2.0?"s":"",
			prefs.psf_step?"manual":"automatic");
*/

/*-- Init the PSF */
    NFPRINTF(OUTPUT,"Initializing PSF modules...");
    psf[ext] = psf_init(prefs.context_name, prefs.context_group,
		prefs.ncontext_name,
		prefs.group_deg, prefs.ngroup_deg,
		set->retisize[0], set->retisize[1],
		psfstep,
		set->nsample);

    psf[ext]->samples_loaded = set->nsample;
    psf[ext]->fwhm = set->fwhm;

/*-- Make the basic PSF-model (1st pass) */
    NFPRINTF(OUTPUT,"Modeling the PSF.");
    psf_make(psf[ext], set);

/*-- Remove bad PSF candidates */
    if (set->nsample>1)
      {
      psf_clean(psf[ext], set, 1);

/*---- Make the basic PSF-model (2nd pass) */
      NFPRINTF(OUTPUT,"Modeling the PSF.");
      psf_make(psf[ext], set);
      }

/*-- Remove bad PSF candidates */
    if (set->nsample>1)
      psf_clean(psf[ext], set, 1);

    psf[ext]->samples_accepted = set->nsample;

/*-- Refine the PSF-model */
    psf_refine(psf[ext], set, prefs.nsuper);

/*-- Just check the Chi2 */
    psf[ext]->chi2 = set->nsample? psf_clean(psf[ext], set, 0) : 0.0;
    NFPRINTF(OUTPUT, "");

/*-- Make a diagnostic of the PSF */
    psf_diagnostic(psf[ext]);
    nmed = ((prefs.context_nsnap-1)/2)*(prefs.context_nsnap+1);
    QPRINTF(OUTPUT, "[%3d/%-3d]     %5d/%-5d  %6.2f    %6.2f %6.2f       %5.3f"
	"      %5.2f\n",
	ext+1, next,
	psf[ext]->samples_accepted, psf[ext]->samples_loaded,
	psfstep,
	psf[ext]->chi2,
	sqrt(psf[ext]->moffat[nmed].fwhm_min*psf[ext]->moffat[nmed].fwhm_max),
	psf[ext]->moffat[nmed].fwhm_max/psf[ext]->moffat[nmed].fwhm_min,
	psf[ext]->moffat[nmed].residuals);

/*-- Load the PCs */
    if (prefs.pc_flag && set->nsample)
      {
      NFPRINTF(OUTPUT,"Including principal components...");
      pc = pc_load(prefs.pc_name);
      pcc = pc_convolve(pc, psf[ext]);
      pco = pc_orthogon(pc, pcc, psf[ext]->pixstep);
      }
    else
      pc = pcc = pco = NULL;

/*-- Save result */
    NFPRINTF(OUTPUT,"Saving the PSF description...");

    psf_save(psf[ext], pco, pc, prefs.psf_name, ext, next);

/*-- Save "Check-images" */
    for (i=0; i<prefs.ncheck_type; i++)
      if (prefs.check_type[i])
        {
        sprintf(str, "Saving CHECK-image #%d...", i+1);
        NFPRINTF(OUTPUT, str);
        psf_writecheck(psf[ext], pco, set, prefs.check_name[i],
		prefs.check_type[i], ext, next);
        }

/*-- Update XML */
    if (prefs.xml_flag)
      update_xml(psf[ext], ncat);

/*-- Free memory */
    end_set(set);

    if (pc)
      {
      pc_end(pc);
      pc_end(pcc);
      pc_end(pco);
      }
    }
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


