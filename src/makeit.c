/*
*				makeit.c
*
* Main loop.
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
#include        "config.h"
#endif

#ifdef USE_THREADS
 #ifdef HAVE_MKL
  #include MKL_H
 #endif
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
#include	"context.h"
#include	"cplot.h"
#include	"diagnostic.h"
#include	"field.h"
#include	"homo.h"
#include	"catout.h"
#include	"pca.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"xml.h"

void		write_error(const char *msg1, const char *msg2);
time_t		thetime, thetime2;

/********************************** makeit ***********************************/
/*
*/
void	makeit(void)

  {
   wcsstruct		*wcs;
   fieldstruct		**fields,
			*field;
   psfstruct		**cpsf,
			*psf;
   setstruct		*set, *set2;
   contextstruct	*context, *fullcontext;
   outcatstruct		*outcat;
   struct tm		*tm;
   char			str[MAXCHAR];
   char			**incatnames,
			*pstr;
   float		**psfbasiss,
			*psfsteps, *psfbasis, *basis,
			psfstep, step;
   int			c,i,p, ncat, ext, next, nmed, nbasis;

/* Install error logging */
  error_installfunc(write_error);

  incatnames = prefs.incat_name;
  ncat = prefs.ncat;

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


/* End here if no filename has been provided */
  if (!ncat)
    {
/*-- Processing end date and time */
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);

/*-- Write XML */
    if (prefs.xml_flag)
      {
      init_xml(0);
      write_xml(prefs.xml_name);
      end_xml();
      }
    return;
    }

/* Create an array of PSFs (one PSF for each extension) */
  QMALLOC(fields, fieldstruct *, ncat);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT, "----- %d input catalogues:\n", ncat);
  for (c=0; c<ncat; c++)
    {
    fields[c] = field_init(incatnames[c]);
    QPRINTF(OUTPUT, "%-20.20s:  \"%-16.16s\"  %3d extension%s %7d detection%s\n",
        fields[c]->rcatname, fields[c]->ident,
        fields[c]->next, fields[c]->next>1 ? "s":"",
        fields[c]->ndet, fields[c]->ndet>1 ? "s":"");
    }
  QPRINTF(OUTPUT, "\n");

  next = fields[0]->next;

  if (prefs.xml_flag)
    init_xml(ncat);  

  psfstep = prefs.psf_step;
  psfsteps = NULL;
  nbasis = 0;
  psfbasis = NULL;
  psfbasiss = NULL;

/* Initialize context */
  NFPRINTF(OUTPUT, "Initializing contexts...");
  context = context_init(prefs.context_name, prefs.context_group,
		prefs.ncontext_group, prefs.group_deg, prefs.ngroup_deg,
		CONTEXT_REMOVEHIDDEN);
  fullcontext = context->npc?
		  context_init(prefs.context_name, prefs.context_group,
		    prefs.ncontext_group, prefs.group_deg, prefs.ngroup_deg,
		    CONTEXT_KEEPHIDDEN)
		: context;

  if (context->npc && ncat<2)
    warning("Hidden dependencies cannot be derived from",
	" a single catalog");
  else if (context->npc && prefs.stability_type == STABILITY_EXPOSURE)
    warning("Hidden dependencies have no effect",
	" in STABILITY_TYPE EXPOSURE mode");

/* Compute PSF steps */
  if (!prefs.psf_step)
    {
    NFPRINTF(OUTPUT, "Computing optimum PSF sampling steps...");
    if (prefs.newbasis_type==NEWBASIS_PCACOMMON
	|| (prefs.stability_type == STABILITY_SEQUENCE
		&& prefs.psf_mef_type == PSF_MEF_COMMON))
      {
      set = load_samples(incatnames, 0, ncat, ALL_EXTENSIONS, next, context);
      psfstep = (float)((set->fwhm/2.35)*0.5);
      end_set(set);
      }
/*-- Need to derive a common pixel step for each ext */
    else if (prefs.newbasis_type == NEWBASIS_PCAINDEPENDENT
	|| context->npc
	|| (prefs.stability_type == STABILITY_SEQUENCE
		&& prefs.psf_mef_type == PSF_MEF_INDEPENDENT))
      {
/*-- Run through all samples to derive a different pixel step for each extension */
      QMALLOC(psfsteps, float, next);
      for (ext=0 ; ext<next; ext++)
        {
        set = load_samples(incatnames, 0, ncat, ext, next, context);
        psfsteps[ext] = (float)(psfstep? psfstep : (set->fwhm/2.35)*0.5);
        end_set(set);
        }
      }
    }

/* Derive a new common PCA basis for all extensions */
  if (prefs.newbasis_type==NEWBASIS_PCACOMMON)
    {
    QMALLOC(cpsf, psfstruct *, ncat*next);
    for (ext=0 ; ext<next; ext++)
      for (c=0; c<ncat; c++)
        {
        sprintf(str, "Computing new PCA image basis from %s...",
		fields[c]->rtcatname);
        NFPRINTF(OUTPUT, str);
        set = load_samples(incatnames, c, 1, ext, next, context);
        step = psfstep;
        cpsf[c+ext*ncat] = make_psf(set, psfstep, NULL, 0, context);
        end_set(set);
        }
    nbasis = prefs.newbasis_number;
    psfbasis = pca_onsnaps(cpsf, ncat*next, nbasis);
    for (i=0 ; i<ncat*next; i++)
      psf_end(cpsf[i]);
    free(cpsf);
    }
/* Derive a new PCA basis for each extension */
  else if (prefs.newbasis_type == NEWBASIS_PCAINDEPENDENT)
    {
    nbasis = prefs.newbasis_number;
    QMALLOC(psfbasiss, float *, next);
    for (ext=0; ext<next; ext++)
      {
      if (psfsteps)
        step = psfsteps[ext];
      else
        step = psfstep;
      QMALLOC(cpsf, psfstruct *, ncat);
      for (c=0; c<ncat; c++)
        {
        sprintf(str, "Computing new PCA image basis from %s...",
		fields[c]->rtcatname);
        NFPRINTF(OUTPUT, str);
        set = load_samples(incatnames, c, 1, ext, next, context);
        cpsf[c] = make_psf(set, step, NULL, 0, context);
        end_set(set);
        }
      psfbasiss[ext] = pca_onsnaps(cpsf, ncat, nbasis);
      for (c=0 ; c<ncat; c++)
        psf_end(cpsf[c]);
      free(cpsf);
      }
    }

  if (context->npc && prefs.hidden_mef_type == HIDDEN_MEF_COMMON)
/*-- Derive principal components of PSF variation from the whole mosaic */
    {
    p = 0;
    QMALLOC(cpsf, psfstruct *, ncat*next);
    for (c=0; c<ncat; c++)
      {
      sprintf(str, "Computing hidden dependency parameter(s) from %s...",
		fields[c]->rtcatname);
      NFPRINTF(OUTPUT, str);
      for (ext=0 ; ext<next; ext++)
        {
        set = load_samples(incatnames, c, 1, ext, next, context);
        if (psfsteps)
          step = psfsteps[ext];
        else
          step = psfstep;
        basis = psfbasiss? psfbasiss[ext] : psfbasis;
        cpsf[p++] = make_psf(set, step, basis, nbasis, context);
        end_set(set);
        }
      }
    free(fullcontext->pc);
    fullcontext->pc = pca_oncomps(cpsf, next, ncat, context->npc);
    for (c=0 ; c<ncat*next; c++)
      psf_end(cpsf[c]);
    free(cpsf);
    }

/* Compute "final" PSF models */
  if (prefs.psf_mef_type == PSF_MEF_COMMON)
    {
    if (prefs.stability_type == STABILITY_SEQUENCE)
      {
/*---- Load all the samples at once */
      set = load_samples(incatnames, 0, ncat, ALL_EXTENSIONS, next, context);
      step = psfstep;
      basis = psfbasis;
      field_count(fields, set, COUNT_LOADED);
      psf = make_psf(set, step, basis, nbasis, context);
      field_count(fields, set, COUNT_ACCEPTED);
      end_set(set);
      NFPRINTF(OUTPUT, "Computing final PSF model...");
      context_apply(context, psf, fields, ALL_EXTENSIONS, 0, ncat);
      psf_end(psf);
      }
    else
      for (c=0; c<ncat; c++)
        {
/*------ Load the samples for current exposure */
        sprintf(str, "Computing final PSF model from %s...",
		fields[c]->rtcatname);
        NFPRINTF(OUTPUT, str);
        set = load_samples(incatnames, c, 1, ALL_EXTENSIONS, next, context);
        if (psfstep)
          step = psfstep;
        else
          step = (float)((set->fwhm/2.35)*0.5);
        basis = psfbasis;
        field_count(fields, set, COUNT_LOADED);
        psf = make_psf(set, step, basis, nbasis, fullcontext);
        field_count(fields, set, COUNT_ACCEPTED);
        end_set(set);
        context_apply(fullcontext, psf, fields, ALL_EXTENSIONS, c, 1);
        psf_end(psf);
        }
    }
  else
    for (ext=0 ; ext<next; ext++)
      {
      basis = psfbasiss? psfbasiss[ext] : psfbasis;
      if (context->npc && prefs.hidden_mef_type == HIDDEN_MEF_INDEPENDENT)
/*------ Derive principal components of PSF components */
        {
        QMALLOC(cpsf, psfstruct *, ncat);
        if (psfsteps)
          step = psfsteps[ext];
        else
          step = psfstep;
        for (c=0; c<ncat; c++)
          {
          if (next>1)
            sprintf(str,
		"Computing hidden dependency parameter(s) from %s[%d/%d]...",
		fields[c]->rtcatname, ext+1, next);
          else
            sprintf(str, "Computing hidden dependency parameter(s) from %s...",
		fields[c]->rtcatname);
          NFPRINTF(OUTPUT, str);
          set = load_samples(incatnames, c, 1, ext, next, context);
          field_count(fields, set, COUNT_LOADED);
          cpsf[c] = make_psf(set, step, basis, nbasis, context);
          field_count(fields, set, COUNT_ACCEPTED);
          end_set(set);
          }
        free(fullcontext->pc);
        fullcontext->pc = pca_oncomps(cpsf, 1, ncat, context->npc);
        for (c=0 ; c<ncat; c++)
          psf_end(cpsf[c]);
        free(cpsf);
        }

      if (prefs.stability_type == STABILITY_SEQUENCE)
        {
/*------ Load all the samples at once */
        if (next>1)
          {
          sprintf(str, "Computing final PSF model for extension [%d/%d]...",
		ext+1, next);
          NFPRINTF(OUTPUT, str);
          }
        else
          NFPRINTF(OUTPUT, "Computing final PSF model...");
        set = load_samples(incatnames, 0, ncat, ext, next, fullcontext);
        if (psfstep)
          step = psfstep;
        else if (psfsteps)
          step = psfsteps[ext];
        else
          step = (float)((set->fwhm/2.35)*0.5);
        field_count(fields, set, COUNT_LOADED);
        psf = make_psf(set, step, basis, nbasis, fullcontext);
        field_count(fields, set, COUNT_ACCEPTED);
        end_set(set);
        context_apply(fullcontext, psf, fields, ext, 0, ncat);
        psf_end(psf);
        }
      else
        for (c=0; c<ncat; c++)
          {
/*-------- Load the samples for current exposure */
          if (next>1)
            sprintf(str, "Reading data from %s[%d/%d]...",
		fields[c]->rtcatname, ext+1, next);
          else
            sprintf(str, "Reading data from %s...",
		fields[c]->rtcatname);
          NFPRINTF(OUTPUT, str);
          set = load_samples(incatnames, c, 1, ext, next, context);
          if (psfstep)
            step = psfstep;
          else if (psfsteps)
            step = psfsteps[ext];
          else
            step = (float)((set->fwhm/2.35)*0.5);
          if (next>1)
            sprintf(str, "Computing final PSF model for %s[%d/%d]...",
		fields[c]->rtcatname, ext+1, next);
          else
            sprintf(str, "Computing final PSF model for %s...",
		fields[c]->rtcatname);
          NFPRINTF(OUTPUT, str);
          field_count(fields, set, COUNT_LOADED);
          psf = make_psf(set, step, basis, nbasis, context);
          field_count(fields, set, COUNT_ACCEPTED);
          end_set(set);
          context_apply(context, psf, fields, ext, c, 1);
          psf_end(psf);
          }
      }

  free(psfsteps);
  if (psfbasiss)
    {
    for (ext=0; ext<next; ext++)
      free(psfbasiss[ext]);
    free(psfbasiss);
    }
  else if (psfbasis)
    free(psfbasis);

/* Compute diagnostics and check-images */
#ifdef USE_THREADS
/* Force MKL using a single thread as diagnostic code is already multithreaded*/ 
 #ifdef HAVE_MKL
  if (prefs.context_nsnap>2)
    mkl_set_num_threads(1);
 #endif
#endif
  QIPRINTF(OUTPUT,
        "   filename      [ext] accepted/total samp. chi2/dof FWHM ellip."
	" resi. asym.");

/* Initialize output catalog */
  if (prefs.outcat_type != CAT_NONE)
    outcat = init_outcat(prefs.outcat_name, context->ncontext);

  for (c=0; c<ncat; c++)
    {
    field = fields[c];
    for (ext=0 ; ext<next; ext++)
      {
      psf = field->psf[ext];
      wcs = field->wcs[ext];
      if (next>1)
        sprintf(str, "Computing diagnostics for %s[%d/%d]...",
		field->rtcatname, ext+1, next);
      else
        sprintf(str, "Computing diagnostics for %s...",
		field->rtcatname);
      NFPRINTF(OUTPUT, str);
/*---- Check PSF with individual datasets */
      set2 = load_samples(incatnames, c, 1, ext, next, context);
      psf->samples_loaded = set2->nsample;
      if (set2->nsample>1)
        {
/*------ Remove bad PSF candidates */
        psf_clean(psf, set2, prefs.prof_accuracy);
        psf->chi2 = set2->nsample? psf_chi2(psf, set2) : 0.0;
        }
      psf->samples_accepted = set2->nsample;
/*---- Compute diagnostics and field statistics */
      psf_diagnostic(psf);
      psf_wcsdiagnostic(psf, wcs);
      nmed = psf->nmed;
      field_stats(fields, set2);
/*---- Display stats for current catalog/extension */
      if (next>1)
        sprintf(str, "[%d/%d]", ext+1, next);
      else
        str[0] = '\0';
      QPRINTF(OUTPUT, "%-17.17s%-7.7s %5d/%-5d %6.2f %6.2f %6.2f  %4.2f"
	" %5.2f %5.2f\n",
	ext==0? field->rtcatname : "",
	str,
	psf->samples_accepted, psf->samples_loaded,
	psf->pixstep,
	psf->chi2,
	psf->moffat_fwhm,
	psf->moffat_ellipticity,
	psf->pfmoffat_residuals,
	psf->sym_residuals);
/*---- Save "Check-images" */
      for (i=0; i<prefs.ncheck_type; i++)
        if (prefs.check_type[i])
          {
          sprintf(str, "Saving CHECK-image #%d...", i+1);
          NFPRINTF(OUTPUT, str);
          check_write(field, set2, prefs.check_name[i], prefs.check_type[i],
		ext, next, prefs.check_cubeflag);
          }
/*---- Write Catalog */
      if (prefs.outcat_type != CAT_NONE)
        write_outcat(outcat, set2);
/*---- Free memory */
      end_set(set2);
      }
    }

  if (prefs.outcat_type != CAT_NONE)
    end_outcat(outcat);

#ifdef USE_THREADS
/* Back to multithreaded MKL */
 #ifdef HAVE_MKL
   mkl_set_num_threads(prefs.nthreads);
 #endif
#endif

/* Save result */
  for (c=0; c<ncat; c++)
    {
    sprintf(str, "Saving PSF model and metadata for %s...",
	fields[c]->rtcatname);
    NFPRINTF(OUTPUT, str);
/*-- Create a file name with a "PSF" extension */
    if (*prefs.psf_dir)
      {
      if ((pstr = strrchr(incatnames[c], '/')))
        pstr++;
      else
        pstr = incatnames[c];
      sprintf(str, "%s/%s", prefs.psf_dir, pstr);
      }
    else
      strcpy(str, incatnames[c]);
    if (!(pstr = strrchr(str, '.')))
      pstr = str+strlen(str);
    sprintf(pstr, "%s", prefs.psf_suffix);
    field_psfsave(fields[c], str);
/* Create homogenisation kernels */
    if (prefs.homobasis_type != HOMOBASIS_NONE)
      {
      for (ext=0; ext<next; ext++)
        {
        if (next>1)
          sprintf(str, "Computing homogenisation kernel for %s[%d/%d]...",
		fields[c]->rtcatname, ext+1, next);
        else
          sprintf(str, "Computing homogenisation kernel for %s...",
		fields[c]->rtcatname);
        NFPRINTF(OUTPUT, str);
        if (*prefs.homokernel_dir)
          {
          if ((pstr = strrchr(incatnames[c], '/')))
            pstr++;
          else
            pstr = incatnames[c];
          sprintf(str, "%s/%s", prefs.homokernel_dir, pstr);
          }
        else
          strcpy(str, incatnames[c]);
        if (!(pstr = strrchr(str, '.')))
          pstr = str+strlen(str);
        sprintf(pstr, "%s", prefs.homokernel_suffix);
        psf_homo(fields[c]->psf[ext], str, prefs.homopsf_params,
		prefs.homobasis_number, prefs.homobasis_scale, ext, next);
        }
      }
#ifdef HAVE_PLPLOT
/* Plot diagnostic maps for all catalogs */
    cplot_ellipticity(fields[c]);
    cplot_fwhm(fields[c]);
    cplot_moffatresi(fields[c]);
    cplot_asymresi(fields[c]);
    cplot_counts(fields[c]);
    cplot_countfrac(fields[c]);
    cplot_modchi2(fields[c]);
    cplot_modresi(fields[c]);
#endif
/*-- Update XML */
    if (prefs.xml_flag)
      update_xml(fields[c]);
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
  for (c=0; c<ncat; c++)
    field_end(fields[c]);
  free(fields);

  if (context->npc)
    context_end(fullcontext);   
  context_end(context);   

  return;
  }


/****** make_psf *************************************************************
PROTO	psfstruct *make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context)
PURPOSE	Make PSFs from a set of FITS binary catalogs.
INPUT	Pointer to a sample set,
	PSF sampling step,
	Pointer to basis image vectors,
	Number of basis vectors,
	Pointer to context structure.
OUTPUT  Pointer to the PSF structure.
NOTES   Diagnostics are computed only if diagflag != 0.
AUTHOR  E. Bertin (IAP)
VERSION 03/09/2009
 ***/
psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context)
  {
   psfstruct		*psf;
   basistypenum		basistype;
   float		pixsize[2];

  pixsize[0] = (float)prefs.psf_pixsize[0];
  pixsize[1] = (float)prefs.psf_pixsize[1];
//  NFPRINTF(OUTPUT,"Initializing PSF modules...");
  psf = psf_init(context, prefs.psf_size, psfstep, pixsize, set->nsample);

  psf->samples_loaded = set->nsample;
  psf->fwhm = set->fwhm;
  
/* Make the basic PSF-model (1st pass) */
//  NFPRINTF(OUTPUT,"Modeling the PSF (1/3)...");
  psf_make(psf, set, 0.2);
  if (basis && nbasis)
    {
    QMEMCPY(basis, psf->basis, float, nbasis*psf->size[0]*psf->size[1]);
    psf->nbasis = nbasis;
    }
  else
    {
//    NFPRINTF(OUTPUT,"Generating the PSF model...");
    basistype = prefs.basis_type;
    if (basistype==BASIS_PIXEL_AUTO)
      basistype = (psf->fwhm < PSF_AUTO_FWHM)? BASIS_PIXEL : BASIS_NONE;
    psf_makebasis(psf, set, basistype, prefs.basis_number);
    }
  psf_refine(psf, set);

/* Remove bad PSF candidates */
  if (set->nsample>1)
    {
    psf_clean(psf, set, 0.2);
 
/*-- Make the basic PSF-model (2nd pass) */
//    NFPRINTF(OUTPUT,"Modeling the PSF (2/3)...");
    psf_make(psf, set, 0.1);
    psf_refine(psf, set);
    }
 
/* Remove bad PSF candidates */
  if (set->nsample>1)
    {
    psf_clean(psf, set, 0.1);
 
/*-- Make the basic PSF-model (3rd pass) */
//    NFPRINTF(OUTPUT,"Modeling the PSF (3/3)...");
    psf_make(psf, set, 0.05);
    psf_refine(psf, set);
    }
 
/* Remove bad PSF candidates */
  if (set->nsample>1)
    psf_clean(psf, set, 0.05);
 
  psf->samples_accepted = set->nsample;
 
/* Refine the PSF-model */
  psf_make(psf, set, prefs.prof_accuracy);
  psf_refine(psf, set);
 
/* Clip the PSF-model */
  psf_clip(psf);

/*-- Just check the Chi2 */
  psf->chi2 = set->nsample? psf_chi2(psf, set) : 0.0;

  return psf;
  }


/****** write_error ********************************************************
PROTO	void    write_error(const char *msg1, const char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/12/2013
 ***/
void	write_error(const char *msg1, const char *msg2)
  {
   char	error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    write_xmlerror(prefs.xml_name, error);
  end_xml();

  return;
  }


