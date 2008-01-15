 /*
				sample.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Read and filter input samples from catalogs.
*
*	Last modify:	15/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "sample.h"
#include "vignet.h"

/******************************** load_samples *******************************/
/*
Examine and load PSF candidates.
*/
setstruct *load_samples(char **filename, int ncat, int ext, int next)
  {
   setstruct		*set;
   catstruct		*cat;
   tabstruct		*tab;
   keystruct		*fkey, *(key[4]);
   char			keynames[4][16]={"FLUX_RADIUS", "FLUX_MAX", "FLAGS",
					"ELONGATION"};
   char			str[MAXCHAR];
   char			*head, *(pkeynames[4]);
   float		*fwhm,*fwhmt,*fwhmt2, *hl, *fmax, *elong,
			backnoise, df,dfmin,fmin, minsn, maxelong, fwhmmin,
			fwhmmax, fval;
   short		*flags;
   int			i,j,n, nobj,nobjmax, imin, nw, ldflag, ext2;

  NFPRINTF(OUTPUT,"Loading samples...");
  minsn = (float)prefs.minsn;
  maxelong = (float)prefs.maxelong;
  fwhmmin = prefs.fwhmrange[0];
  fwhmmax = prefs.fwhmrange[1];
  fwhm = NULL;	/* To avoid gcc -Wall warnings */

  if (prefs.autoselect_flag)
    {
/*-- Allocate memory */
    nobj = 0;
    nobjmax = LSAMPLE_DEFSIZE;
    QMALLOC(fwhm, float, nobjmax);
    fwhmt=fwhm;

/*-- Initialize string array */
    for (i=0; i<4; i++)
      pkeynames[i] = keynames[i];

/*-- Try to estimate the most appropriate Half-light Radius range */
/*-- Get the Half-light radii */
    nobj = 0;
    for (i=0; i<ncat; i++)
      {
      sprintf(str,"Examining Catalog #%d", i+1);
      NFPRINTF(OUTPUT, str);
/*---- Read input catalog */
      if (!(cat = read_cat(filename[i])))
        error(EXIT_FAILURE, "*Error*: No such catalog: ", filename[i]);

/*---- Get the background noise for this catalog */
      ldflag = 1;
      ext2 = ext+1;
      tab = cat->tab;
      for (j=cat->ntab; j--; tab=tab->nexttab)
        if (!(ldflag = strcmp("LDAC_IMHEAD",tab->extname))
	  || fitsread(tab->headbuf,"SEXBKDEV",&backnoise,H_FLOAT,T_FLOAT)
		==RETURN_OK)
          if (!--ext2)
            break;

      if (j<0)
        error(EXIT_FAILURE, "*Error*: SExtractor table missing in ",
		filename[j]);
      if (!ldflag)
        {
        fkey=read_key(tab, "Field Header Card");
        head = fkey->ptr;
        if (fitsread(head,"SEXBKDEV",&backnoise,H_FLOAT,T_FLOAT)==RETURN_ERROR)
          error(EXIT_FAILURE, "*Error*: Keyword not found:", "SEXBKDEV");
	}
      if (backnoise<1/BIG)
        backnoise = 1.0;

/*---- Now load the objects */
      ext2 = ext+1;
      tab = cat->tab;
      for (j=cat->ntab; j--; tab=tab->nexttab)
        if (!strcmp("LDAC_OBJECTS", tab->extname)
		||  !strcmp("OBJECTS", tab->extname))
          if (!--ext2)
            break;
      if (j<0)
        error(EXIT_FAILURE, "*Error*: OBJECTS table not found in catalog ",
		filename[i]);
      for (j=0; j<4; j++)
        if (!(key[j]=name_to_key(tab, keynames[j])))
          {
          sprintf(str, "%s not found in catalog %s", keynames[j], filename[i]);
          error(EXIT_FAILURE, "*Error*: ", str);
          }

      read_keys(tab, pkeynames, NULL, 4, NULL);

/*---- Fill the FWHM array */
      hl = key[0]->ptr;
      fmax = key[1]->ptr;
      flags = key[2]->ptr;
      elong = key[3]->ptr;
      for (n=tab->naxisn[1]; n--; hl++, fmax++, flags++, elong++)
	{
        if (*fmax/backnoise>minsn
		&& !(*flags&prefs.flag_mask)
		&& *elong<maxelong
		&& (fval=2.0**hl)>=fwhmmin
		&& fval<fwhmmax)
          {
          if (++nobj>nobjmax)
            {
            nobjmax += LSAMPLE_DEFSIZE;
            QREALLOC(fwhm, float, nobjmax);
            fwhmt=fwhm+nobj-1;
            }
          *(fwhmt++) = fval;
          }
	}
      free_cat(&cat, 1);
      }

    if (nobj)
      {
/*---- Sort FWHMs */
      hmedian(fwhm, nobj);

/*---- Find the mode */
      nw = nobj/4;
      if (nw<4)
        nw = 1;
      dfmin = BIG;
      fmin = 0.0;
      imin = 0;
      fwhmt = fwhm;
      fwhmt2 = fwhm+nw;
      for (i=nobj-nw; i--; fwhmt++,fwhmt2++)
        {  
        if ((df = *fwhmt2 - *fwhmt) < dfmin)
          {
          dfmin = df;
          fmin = (*fwhmt2 + *fwhmt)/2.0;
          imin++;
          }
        }

      if (nobj<2)
        fmin = *fwhm;

      dfmin = (float)pow((double)prefs.maxvar+1.0, 0.3333333);
      fwhmmin = dfmin>0.0?fmin/dfmin:0.0;
      if (fwhmmin<prefs.fwhmrange[0])
        fwhmmin = prefs.fwhmrange[0];
      fwhmmax = fmin*dfmin*dfmin;
      if (fwhmmax>prefs.fwhmrange[1])
        fwhmmax = prefs.fwhmrange[1];
      }
    else
      fmin = 2.35/(1.0-1.0/INTERPFAC);
    }
  else
    {
    warning("No source with appropriate FWHM found!!","");
    fmin = fwhmmin = fwhmmax = 0.0;
    }

  free(fwhm);
/*
  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "Selected FWHM range: %.2f - %.2f pixels\n",
		fwhmmin, fwhmmax);
*/
/* Load the samples */
  set = NULL;
  for (i=0; i<ncat; i++)
    set = read_samples(set, filename[i], fwhmmin/2.0, fwhmmax/2.0, ext, next);

  set->fwhm = fmin;
  sprintf(str, "%d samples loaded.", set->nsample);
  NFPRINTF(OUTPUT, str);
  if (!set->nsample)
    warning("No appropriate source found!!","");

  return set;
  }


/******************************** read_samples *******************************/
/*
*/
setstruct *read_samples(setstruct *set, char *filename,
			float frmin, float frmax, int ext, int next)

  {
   catstruct		*cat;
   tabstruct		*tab, *keytab;
   keystruct		*key, *vigkey;
   samplestruct		*sample;
   t_type		contexttyp[MAXCONTEXT];
   void			*context[MAXCONTEXT];
   static char		str[MAXCHAR], str2[MAXCHAR];
   char			**kstr,
			*head, *buf;
   unsigned short	*flags;
   double		contextval[MAXCONTEXT], *cmin, *cmax, dval, sn;
   float		*xm, *ym, *vignet,*vignett, *flux, *fluxmax, *fluxrad,
			*elong,
			backnoise, backnoise2, gain, minsn,maxelong;
   static int		ncat;
   int			i,j, n, nsample,nsamplemax,
			vigw, vigh, vigsize, imaw,imah, nobj, nt,
			maxbad, maxbadflag, ldflag, ext2;


  maxbad = prefs.badpix_nmax;
  maxbadflag = prefs.badpix_flag;
  maxelong = prefs.maxelong;
  minsn = prefs.minsn;

/* If a NULL pointer is provided, we allocate a new set */
  if (!set)
    {
    set = init_set();
    nsample = nsamplemax = 0;
    ncat = 1;
    }
  else
    nsample = nsamplemax = set->nsample;

  cmin = cmax = (double *)NULL;	/* To avoid gcc -Wall warnings */
  if (set->ncontext)
    {
    QMALLOC(cmin, double, set->ncontext);
    QMALLOC(cmax, double, set->ncontext);
    for (i=0; i<set->ncontext; i++)
      if (ncat>1 && set->nsample)
        {
        cmin[i] = set->contextoffset[i] - set->contextscale[i]/2.0;
        cmax[i] = cmin[i] + set->contextscale[i];
        }
      else
        {
        cmin[i] = BIG;
        cmax[i] = -BIG;
        }
    }

/*-- Read input catalog */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: No such catalog: ", filename);
  head = (char *)NULL;	/* To avoid gcc -Wall warnings */
  ldflag = 1;
  ext2 = ext+1;
  tab = cat->tab;
  for (j=cat->ntab; j--; tab=tab->nexttab)
    if (!(ldflag = strcmp("LDAC_IMHEAD",tab->extname))
	  || fitsread(head=tab->headbuf,"SEXBKDEV",&backnoise,H_FLOAT,T_FLOAT)
		==RETURN_OK)
      if (!--ext2)
        break;
  if (j<0)
    error(EXIT_FAILURE, "*Error*: SExtractor table missing in ", filename);
  if (!ldflag)
    {
    key=read_key(tab, "Field Header Card");
    head = key->ptr;
    if (fitsread(head,"SEXBKDEV",&backnoise,H_FLOAT,T_FLOAT)==RETURN_ERROR)
      error(EXIT_FAILURE, "*Error*: Keyword not found:", "SEXBKDEV");
    }
  backnoise2 = backnoise*backnoise;
  if ((n=fitsfind(head, "END     ")) != RETURN_ERROR)
    {
    QCALLOC(set->head, char, ((n*80)/FBSIZE+1)*FBSIZE);
    memcpy(set->head, head, (n+1)*80);
    }
  if (fitsread(head, "SEXGAIN", &gain, H_FLOAT, T_FLOAT) == RETURN_ERROR)
    error(EXIT_FAILURE, "*Error*: Keyword not found:", "SEXGAIN");

  if (fitsread(head,  "SEXIMASX", &imaw, H_INT, T_LONG)== RETURN_ERROR
	&& fitsread(head, "NAXIS1  ", &imaw, H_INT, T_LONG)== RETURN_ERROR)
    error(EXIT_FAILURE,"*Error*: Image X-size not found","");

  if (fitsread(head,  "SEXIMASY", &imah, H_INT, T_LONG)== RETURN_ERROR
	&& fitsread(head, "NAXIS2  ", &imah, H_INT, T_LONG)== RETURN_ERROR)
    error(EXIT_FAILURE,"*Error*: Image Y-size not found","");

  ext2 = ext+1;
  tab = cat->tab;
  for (j=cat->ntab; j--; tab=tab->nexttab)
    if (!strcmp("LDAC_OBJECTS", tab->extname)
		||  !strcmp("OBJECTS", tab->extname))
      if (!--ext2)
        break;
  if (j<0)
    error(EXIT_FAILURE, "*Error*: OBJECTS table not found in catalog ",
		filename);

/* Init the single-row tab */
  keytab = init_readobj(tab, &buf);

  if (!(key = name_to_key(keytab, "X_IMAGE")))
    error(EXIT_FAILURE, "*Error*: X_IMAGE parameter not found in catalog ",
		filename);
  xm = (float *)key->ptr;
  if (!(key = name_to_key(keytab, "Y_IMAGE")))
    error(EXIT_FAILURE, "*Error*: Y_IMAGE parameter not found in catalog ",
		filename);
  ym = (float *)key->ptr;

  if (!(key = name_to_key(keytab, "FLUX_RADIUS")))
    error(EXIT_FAILURE, "*Error*: FLUX_RADIUS parameter not found in catalog ",
		filename);
  fluxrad = (float *)key->ptr;

  if (!(key = name_to_key(keytab, "FLUX_APER")))
    error(EXIT_FAILURE, "*Error*: FLUX_APER parameter not found in catalog ",
		filename);
  flux = (float *)key->ptr;

  if (!(key = name_to_key(keytab, "FLUX_MAX")))
    error(EXIT_FAILURE,"*Error*: FLUX_MAX parameter not found in catalog ",
		filename);
  fluxmax = (float *)key->ptr;

  if (!(key = name_to_key(keytab, "ELONGATION")))
    error(EXIT_FAILURE, "*Error*: ELONGATION parameter not found in catalog ",
		filename);
  elong = (float *)key->ptr;

  if (!(key = name_to_key(keytab, "FLAGS")))
    error(EXIT_FAILURE, "*Error*: FLAGS parameter not found in catalog ",
		filename);
  flags = (unsigned short *)key->ptr;
  nobj = key->nobj;

  if (!(key = name_to_key(keytab, "VIGNET")))
    error(EXIT_FAILURE,
	"*Error*: VIGNET parameter not found in catalog ", filename);
  vignet = (float *)key->ptr;
  if (key->naxis != 2)
    error(EXIT_FAILURE, "*Error*: VIGNET should be a 2D vector", "");
  vigkey = key;
  vigw = *(vigkey->naxisn);
  vigh = *(vigkey->naxisn+1);
  vigsize = vigw*vigh;
  if (!set->sample)
    {
    set->vigsize[0] = vigw;
    set->vigsize[1] = vigh;
    set->nvig = vigw*vigh;
    }

/* Try to load the set of context keys */
  kstr = prefs.context_name;
  for (i=0; i<set->ncontext; i++, kstr++)
    if (**kstr==(char)':')
      {
      context[i] = &contextval[i];
      contexttyp[i] = T_DOUBLE;
      if (fitsread(head, *kstr+1, context[i], H_FLOAT,T_DOUBLE)==RETURN_ERROR)
        {
        sprintf(str, "*Error*: %s parameter not found in the header of ",
		*kstr+1);
        error(EXIT_FAILURE, str, filename);
        }
      }
    else
      {
      if (!(key = name_to_key(keytab, *kstr)))
        {
        sprintf(str, "*Error*: %s parameter not found in catalog ", *kstr);
        error(EXIT_FAILURE, str, filename);
        }
      context[i] = key->ptr;
      contexttyp[i] = key->ttype;
      strcpy(set->contextname[i], key->name);
      }
  if (next>1)
    sprintf(str2, "[%d/%d]", ext+1, next);
  else
    strcpy(str2, "");

/* Now examine each vector of the shipment */
  nt = keytab->naxisn[1];
  for (n=0; nt; n++)
    {
    nt = read_obj(keytab,tab, buf);
    if (!(n%100))
      {
      sprintf(str,"Catalog #%d %s: Object #%d / %d samples stored",
	ncat, str2, n,nsample);
      NFPRINTF(OUTPUT, str);
      }
    sn = (double)(backnoise>0.0? *fluxmax/backnoise : BIG);
/*---- Apply some selection over flags, fluxes... */
    if (!(*flags&prefs.flag_mask)
	&& sn>minsn
	&& *fluxrad>frmin && *fluxrad<frmax
	&& *elong<maxelong)
      {
/*---- ... and check the integrity of the sample */
      j = 0;
      vignett = vignet;
      for (i=vigsize; i--; vignett++)
        if (*vignett <= -BIG)
          j++;
      if (maxbadflag && j > maxbad)
        continue; 
    
/*---- Allocate memory for the first shipment */
      if (!set->sample)
        {
/*------ The retina cannot be larger than the input images! */
/*
        if (prefs.retisize[0]*prefs.psf_step>(double)vigw
	 || prefs.retisize[1]*prefs.psf_step>(double)vigh)
          error(EXIT_FAILURE, "*Error*: Vignets smaller than the retina in ",
		filename);
*/
        nsample = 0;
        nsamplemax = LSAMPLE_DEFSIZE;
        malloc_samples(set, nsamplemax);
        }
      else
        {
        if (set->vigsize[0] != vigw || set->vigsize[1] != vigh)
          error(EXIT_FAILURE, "*Error*: Incompatible VIGNET size found in ",
		filename);
        }

/*---- Increase storage space to receive new candidates if needed */
      if (nsample>=nsamplemax)
        {
         int	nadd=(int)(1.62*nsamplemax);
        nsamplemax = nadd>nsamplemax?nadd:nsamplemax+1;
        realloc_samples(set, nsamplemax);
        }

      sample = set->sample + nsample;

/*---- Copy the vignet to the training set */
      memcpy(sample->vig, vignet, vigsize*sizeof(float));

      sample->norm = *flux;
      sample->backnoise2 = backnoise2;
      sample->gain = gain;
/*---- Use a first approximation of the center for feeding the retina */
      sample->x = *xm;
      sample->y = *ym;
      sample->dx = sample->x - (int)(sample->x+0.49999);
      sample->dy = sample->y - (int)(sample->y+0.49999);
      for (i=0; i<set->ncontext; i++)
        {
        dval = sample->context[i];
        ttypeconv(context[i], &dval, contexttyp[i], T_DOUBLE);
        sample->context[i] = dval;
/*------ Update min and max */
        if (dval<cmin[i])
          cmin[i] = dval;
        if (dval>cmax[i])
          cmax[i] = dval;
        }
      make_weights(set, sample);
      nsample++;
      }
    }

/* Update the scaling */
  if (set->ncontext)
    {
    if (nsample)
      for (i=0; i<set->ncontext; i++)
        {
        set->contextscale[i] = cmax[i] - cmin[i];
        set->contextoffset[i] = (cmin[i] + cmax[i])/2.0;
        }
    free(cmin);
    free(cmax);
    }
  end_readobj(keytab,tab, buf);
  free_cat(&cat, 1); 

  set->nsample = nsample;

/* Don't waste memory! */
  if (nsample)
    realloc_samples(set, nsample);

/* Increase the current catalog number */
  ncat++;

  return set;
  }


/****** malloc_samples *******************************************************
PROTO   void malloc_samples(setstruct *set, int nsample)
PURPOSE Allocate memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 02/03/99
*/
void	malloc_samples(setstruct *set, int nsample)

  {
   samplestruct	*sample;
   int		n;

  QMALLOC(set->sample, samplestruct, nsample);
  sample = set->sample;
  for (n=nsample; n--; sample++)
    {
    QMALLOC(sample->vig, float, set->nvig);
    QMALLOC(sample->vigresi, float, set->nvig);
    QMALLOC(sample->vigweight, float, set->nvig);
    QMALLOC(sample->retina, float, set->nreti);
    QMALLOC(sample->retiweight, float, set->nreti);
    if (set->ncontext)
      QMALLOC(sample->context, double, set->ncontext);
    }

  set->nsamplemax = nsample;

  return;
  }


/****** realloc_samples ******************************************************
PROTO   void realloc_samples(setstruct *set, int nsample)
PURPOSE Re-allocate memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 02/03/99
*/
void	realloc_samples(setstruct *set, int nsample)

  {
   samplestruct	*sample;
   int		n;

/* If we want to reallocate 0 samples, better free the whole thing! */
  if (!nsample)
    free_samples(set);

/* Two cases: either more samples are required, or the opposite! */
  if (nsample>set->nsamplemax)
    {
    QREALLOC(set->sample, samplestruct, nsample);
    sample = set->sample + set->nsamplemax;
    for (n = nsample - set->nsamplemax; n--; sample++)
      {
      QMALLOC(sample->vig, float, set->nvig);
      QMALLOC(sample->vigresi, float, set->nvig);
      QMALLOC(sample->vigweight, float, set->nvig);
      QMALLOC(sample->retina, float, set->nreti);
      QMALLOC(sample->retiweight, float, set->nreti);
      if (set->ncontext)
        QMALLOC(sample->context, double, set->ncontext);
      }
    }
  else if (nsample<set->nsamplemax)
    {
    sample = set->sample + nsample;
    for (n = set->nsamplemax - nsample; n--; sample++)
      {
      free(sample->vig);
      free(sample->vigresi);
      free(sample->vigweight);
      free(sample->retina);
      free(sample->retiweight);
      if (set->ncontext)
        free(sample->context);
      }
    QREALLOC(set->sample, samplestruct, nsample);
    }

  set->nsamplemax = nsample;

  return;
  }


/****** free_samples *********************************************************
PROTO   void free_samples(setstruct *set, int nsample)
PURPOSE free memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 02/03/99
*/
void	free_samples(setstruct *set)

  {
   samplestruct	*sample;
   int		n;

  sample = set->sample;
  for (n = set->nsamplemax; n--; sample++)
    {
    free(sample->vig);
    free(sample->vigresi);
    free(sample->vigweight);
    free(sample->retina);
    free(sample->retiweight);
    if (set->ncontext)
      free(sample->context);
    }

  free(set->sample);
  set->sample = NULL;
  set->nsample = set->nsamplemax = 0;

  return;
  }


/****** remove_sample ********************************************************
PROTO   samplestruct *remove_sample(setstruct *set, int isample)
PURPOSE Remove an element from a set of samples.
INPUT   set structure pointer,
        sample number.
OUTPUT  The new pointer for the element that replaced the removed one.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 01/03/99
*/
samplestruct	*remove_sample(setstruct *set, int isample)

  {
   static samplestruct	exsample;
   samplestruct		*sample;
   int			nsample;

/* If we want to reallocate 0 samples, better free the whole thing! */
  nsample = set->nsample-1;
  if (nsample>0)
    {
    sample = set->sample + isample;
    exsample = *(set->sample+nsample);
    *(set->sample+nsample) = *sample;
    *sample = exsample;
    }
   else
     nsample=0;
  realloc_samples(set, nsample);
  set->nsample = nsample;

  return set->sample+isample;
  }


/****** init_set ************************************************************
PROTO   setstruct *init_set()
PURPOSE Allocate and initialize a set structure.
INPUT   -.
OUTPUT  -.
NOTES   See prefs.h.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 25/04/2007
*/
setstruct	*init_set(void)

  {
   setstruct	*set;
   int		i;

  QCALLOC(set, setstruct, 1);
  set->nsample = set->nsamplemax = 0;
  set->vigdim = set->retidim = 2;
  QMALLOC(set->retisize, int, set->retidim);
  QMALLOC(set->vigsize, int, set->vigdim);
  set->retisize[0] = prefs.retisize[0];
  set->retisize[1] = prefs.retisize[1];
  set->nreti = set->retisize[0]*set->retisize[1];/* Temporary solution (?) */
  set->ncontext = prefs.ncontext_group;
  if (set->ncontext)
    {
    QMALLOC(set->contextoffset, double, set->ncontext);
    QMALLOC(set->contextscale, double, set->ncontext);
    QMALLOC(set->contextname, char *, set->ncontext);
    for (i=0; i<set->ncontext; i++)
      QMALLOC(set->contextname[i], char, 80);
    }

  return set;
  }


/****** end_set *************************************************************
PROTO   void end_set(setstruct *set)
PURPOSE free memory allocated by a complete set structure.
INPUT   set structure pointer,
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 10/07/2003
*/
void	end_set(setstruct *set)

  {
   int	i;

  free_samples(set);
  free(set->vigsize);
  free(set->retisize);
  if (set->ncontext)
    {
    for (i=0; i<set->ncontext; i++)
      free(set->contextname[i]);
    free(set->contextname);
    free(set->contextoffset);
    free(set->contextscale);
    }
  free(set->head);
  free(set);

  return;
  }


/****** make_weights *********************************************************
PROTO   void make_weights(setstruct *set, samplestruct *sample)
PURPOSE Produce a weight-map for each sample vignet.
INPUT   set structure pointer,
        sample structure pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP,Leiden observatory & ESO)
VERSION 13/08/2007
*/
void make_weights(setstruct *set, samplestruct *sample)

  {
   float	*vig, *vigweight,
		backnoise2, gain, noise2, profaccu2, pix;
   int		i;

/* Produce a weight-map */
  profaccu2 = prefs.prof_accuracy*prefs.prof_accuracy;
  gain = sample->gain;
  backnoise2 = sample->backnoise2;
  for (vig=sample->vig, vigweight=sample->vigweight, i=set->nvig; i--;)
    {
    if (*vig <= -BIG)
      *(vig++) = *(vigweight++) = 0.0;
    else
      {
      pix = *(vig++);
      noise2 = backnoise2 + profaccu2*pix*pix;
      if (pix>0.0 && gain>0.0)
        noise2 += pix/gain;
      *(vigweight++) = 1.0/noise2;      
      }
    }

  return;
  }


/****** update_retina ********************************************************
PROTO   void update_retina(setstruct *set, samplestruct *sample,
                           float dx, float dy, float pixstep)
PURPOSE Update the retina content, copying data from the vignet, normalizing
        it and producing a weight-map.
INPUT   set structure pointer,
        sample structure pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP,Leiden observatory & ESO)
VERSION 05/08/99
*/
void update_retina(setstruct *set, samplestruct *sample, float pixstep)

  {
   float	*retina, *retiweight,
		backnoise2, gain, norm, norm2, noise2, profaccu2, pix;
   int		i;

  vignet_resample(sample->vig, set->vigsize[0], set->vigsize[1],
	sample->retina, set->retisize[0], set->retisize[1],
	sample->dx, sample->dy, pixstep, pixstep>1.0?pixstep:1.0);
	  
/* Normalize approximately the retina and produce a weight-map */
  norm = sample->norm;
  norm2 = norm*norm;
  profaccu2 = prefs.prof_accuracy*prefs.prof_accuracy*norm2;
  gain = sample->gain;
  backnoise2 = sample->backnoise2;
  retina = sample->retina;
  retiweight = sample->retiweight;
  for (i=set->nreti; i--;)
    {
    pix = (*(retina++) /= norm);
    noise2 = backnoise2 + profaccu2*pix*pix;
    if (pix>0.0 && gain>0.0)
      noise2 += pix/gain;
    *(retiweight++) = norm2/noise2;      
    }

  return;
  }


