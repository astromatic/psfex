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
*	Last modify:	31/10/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"check.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"

/********************************** makeit ***********************************/
/*
*/
void	makeit(char **incatnames, int ncat)

  {
   setstruct		*set;
   psfstruct		*psf;
   pcstruct		*pc, *pcc, *pco;
   catstruct		*cat;
   tabstruct		*tab;
   static char		str[MAXCHAR];
   float		psfstep;
   int			i, ntab, ext, next;

   out_data_struct  out_data;

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
  for (ext = 0 ; ext<next; ext++)
    {
/*-- Load all the samples */
    NFPRINTF(OUTPUT,"Loading samples...");
    set = load_samples(incatnames, ncat, ext, next);

    NFPRINTF(OUTPUT, "");
    NPRINTF(OUTPUT, "%d samples loaded\n", set->nsample);
    out_data.samples_loaded = set->nsample;

    if (!set->nsample)
      warning("No appropriate source found!!","");

    psfstep = (float)(prefs.psf_step? prefs.psf_step
				: (set->fwhm/2.35)*(1.0-1.0/INTERPFAC));
    out_data.fwhm = set->fwhm;

    NFPRINTF(OUTPUT, "");
    NPRINTF(OUTPUT, "PSF sampled each %.2f pixel%s (%s)\n",
			psfstep,
			psfstep>=2.0?"s":"",
			prefs.psf_step?"manual":"automatic");

/*-- Init the PSF */
    NFPRINTF(OUTPUT,"Initializing PSF modules...");
    psf = psf_init(prefs.context_name, prefs.context_group,prefs.ncontext_name,
		prefs.group_deg, prefs.ngroup_deg,
		set->retisize[0], set->retisize[1],
		psfstep,
		set->nsample);

/*-- Make the basic PSF-model (1st pass) */
    NFPRINTF(OUTPUT,"Modeling the PSF.");
    psf_make(psf, set);

/*-- Remove bad PSF candidates */
    if (set->nsample>1)
      {
      psf_clean(psf, set, 1);
      NFPRINTF(OUTPUT, "");
      NPRINTF(OUTPUT, "%d samples accepted\n", set->nsample);

/*---- Make the basic PSF-model (2nd pass) */
      NFPRINTF(OUTPUT,"Modeling the PSF.");
      psf_make(psf, set);
      }

/*-- Remove bad PSF candidates */
    if (set->nsample>1)
      {
      psf_clean(psf, set, 1);
      NFPRINTF(OUTPUT, "");
      NPRINTF(OUTPUT, "%d samples accepted\n", set->nsample);
      }

    out_data.samples_accepted = set->nsample;

/*-- Refine the PSF-model */
    psf_refine(psf, set, prefs.nsuper);

/*-- Just check the Chi2 */
    out_data.chi2 = set->nsample? psf_clean(psf, set, 0) : 0.0;
    NFPRINTF(OUTPUT, "");

/*-- Load the PCs */
    if (prefs.pc_flag && set->nsample)
      {
      NFPRINTF(OUTPUT,"Including principal components...");
      pc = pc_load(prefs.pc_name);
      pcc = pc_convolve(pc, psf);
      pco = pc_orthogon(pc, pcc, psf->pixstep);
      }
    else
      pc = pcc = pco = NULL;

/*-- Save result */
    NFPRINTF(OUTPUT,"Saving the PSF description...");

    psf_save(psf, pco, pc, prefs.psf_name, ext, next, &out_data);

/*-- Save "Check-images" */
    for (i=0; i<prefs.ncheck_type; i++)
      if (prefs.check_type[i])
        {
        sprintf(str, "Saving CHECK-image #%d...", i+1);
        NFPRINTF(OUTPUT, str);
        psf_writecheck(psf, pco, set, prefs.check_name[i],prefs.check_type[i],
		ext, next);
        }

/*-- Free memory */
    end_set(set);
    psf_end(psf);

    if (pc)
      {
      pc_end(pc);
      pc_end(pcc);
      pc_end(pco);
      }
    }

  return;
  }

