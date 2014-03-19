/*
*				sample.c
*
* Read and filter input samples from catalogues.
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "context.h"
#include "misc.h"
#include "sample.h"
#include "vignet.h"

static float	compute_fwhmrange(float *fwhm, int nfwhm, float maxvar,
		float minin, float maxin, float *minout, float *maxout);

/****** load_samples *********************************************************
PROTO	setstruct *load_samples(char **filename, int catindex, int ncat,
		int ext, int next, contextstruct *context)
PURPOSE	Examine and load point source data.
INPUT	Array of catalog filenames,
	catalog index,
	current extension,
	number of extensions,
	pointer to the context.
OUTPUT  Pointer to a set containing samples that match acceptance criteria.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/07/2012
*/
setstruct *load_samples(char **filename, int catindex, int ncat, int ext,
			int next, contextstruct *context)
  {
   setstruct		*set;
   catstruct		*cat;
   tabstruct		*tab;
   keystruct		*fkey, *key;
   char			keynames[][32]={"ELONGATION", "FLAGS", "FLUX_RADIUS",
				"SNR_WIN", ""};
   char			str[MAXCHAR];
   char			**pkeynames,
			*head;
   float		*fwhmmin,*fwhmmax,*fwhmmode,
			*fwhm,*fwhmt, *elong, *hl, *snr,
			backnoise, minsn, maxelong, min,max, mode,  fval;
   unsigned int		*imaflags;
   unsigned short	*flags, *wflags;
   int			*fwhmindex,
			e,i,j,n, icat, nobj,nobjmax, ldflag, ext2, nkeys;

//  NFPRINTF(OUTPUT,"Loading samples...");
  minsn = (float)prefs.minsn;
  maxelong = (float)(prefs.maxellip < 1.0?
	(prefs.maxellip + 1.0)/(1.0 - prefs.maxellip)
	: 100.0);
  min = (float)prefs.fwhmrange[0];
  max = (float)prefs.fwhmrange[1];
  fwhm = NULL;	/* To avoid gcc -Wall warnings */
/* Allocate memory */
  QMALLOC(fwhmmin, float, ncat);
  QMALLOC(fwhmmax, float, ncat);
  QMALLOC(fwhmmode, float, ncat);

  if (prefs.autoselect_flag)
    {
/*-- Allocate memory */
    nobjmax = LSAMPLE_DEFSIZE;
    QMALLOC(fwhm, float, nobjmax);
    QMALLOC(fwhmindex, int, ncat+1);
    fwhmindex[0] = nobj = 0;
    fwhmt=fwhm;

/*-- Initialize string array */
    for (i=0; (*keynames[i]); i++);
    nkeys = i;
    QMALLOC(pkeynames, char *, nkeys);
    for (i=0; i<nkeys; i++)
      pkeynames[i] = keynames[i];

/*-- Try to estimate the most appropriate Half-light Radius range */
/*-- Get the Half-light radii */
    nobj = 0;
    for (i=0; i<ncat; i++)
      {
      sprintf(str,"Examining Catalog #%d", i+1);
//      NFPRINTF(OUTPUT, str);
/*---- Read input catalog */
      icat = catindex + i;
      if (!(cat = read_cat(filename[icat])))
        error(EXIT_FAILURE, "*Error*: No such catalog: ", filename[icat]);

/*---- Load the objects */
      e = 0;
      e = ext;
      tab = cat->tab;
      for (j=cat->ntab; j--; tab=tab->nexttab)
        if (!strcmp("LDAC_OBJECTS", tab->extname)
		||  !strcmp("OBJECTS", tab->extname))
          {
          if (ext != ALL_EXTENSIONS)
            {
            if (e>0)
              {
              e--;
              continue;
              }
            else if (e<0)
              break;
            e--;
            }

/*-------- Read the data */

          read_keys(tab, pkeynames, NULL, nkeys, NULL);

          if ((key = name_to_key(tab, "ELONGATION")))
            elong = (float *)key->ptr;
          else
            {
            warning("ELONGATION parameter not found in catalog ",
			filename[icat]);
            elong = NULL;
            }
          if ((key = name_to_key(tab, "FLAGS")))
            flags = (unsigned short *)key->ptr;
          else
            {
            warning("FLAGS parameter not found in catalog ", filename[icat]);
            flags = NULL;
            }
          if ((key = name_to_key(tab, "FLAGS_WEIGHT")))
            wflags = (unsigned short *)key->ptr;
          else
            wflags = NULL;
          if ((key = name_to_key(tab, "IMAFLAGS_ISO")))
            imaflags = (unsigned int *)key->ptr;
          else
            imaflags = NULL;
          if (!(key = name_to_key(tab, "FLUX_RADIUS")))
              {
              sprintf(str, "FLUS_RADIUS not found in catalog %s",
			filename[icat]);
              error(EXIT_FAILURE, "*Error*: ", str);
              }
          hl = (float *)key->ptr;
          if (!(key = name_to_key(tab, "SNR_WIN")))
              {
              sprintf(str, "SNR_WIN not found in catalog %s",
			filename[icat]);
              error(EXIT_FAILURE, "*Error*: ", str);
              }
          snr = (float *)key->ptr;

          for (n=tab->naxisn[1]; n--; hl++, snr++, flags++, elong++)
            {
            if (*snr>minsn
		&& !(flags && (*flags&prefs.flag_mask))
		&& !(wflags && (*wflags&prefs.wflag_mask))
		&& !(imaflags && (*imaflags&prefs.imaflag_mask))
		&& (!(elong && *elong>=maxelong))
		&& (fval=2.0**hl)>=min
		&& fval<max)
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
          }
      free_cat(&cat, 1);
      fwhmindex[i+1] = nobj;
      }

    if (prefs.var_type == VAR_NONE)
      {
      if (nobj)
        mode = compute_fwhmrange(fwhm, nobj, prefs.maxvar,
		prefs.fwhmrange[0],prefs.fwhmrange[1], &min, &max);
      else
        {
        warning("No source with appropriate FWHM found!!","");
        mode = min = max = 2.35/(1.0-1.0/INTERPFAC);
        }
      for (i=0; i<ncat; i++)
        {
        fwhmmin[i] = min;
        fwhmmax[i] = max;
        fwhmmode[i] = mode;
        }
      }
    else
      for (i=0; i<ncat; i++)
        {
        nobj = fwhmindex[i+1] - fwhmindex[i];
        if (nobj)
          {
          fwhmmode[i] = compute_fwhmrange(&fwhm[fwhmindex[i]],
		fwhmindex[i+1]-fwhmindex[i], prefs.maxvar,
		prefs.fwhmrange[0],prefs.fwhmrange[1], &fwhmmin[i],&fwhmmax[i]);
          }
        else
          {
          warning("No source with appropriate FWHM found!!","");
          fwhmmode[i] = fwhmmin[i] = fwhmmax[i] = 2.35/(1.0-1.0/INTERPFAC);
          }
        }
    free(fwhm);
    free(fwhmindex);
    free(pkeynames);
    }
  else
    for (i=0; i<ncat; i++)
      {
      fwhmmin[i] = (float)prefs.fwhmrange[0];
      fwhmmax[i] = (float)prefs.fwhmrange[1];
      fwhmmode[i] = (fwhmmin[i] + fwhmmax[i]) / 2.0;
      }


/* Load the samples */
  set = NULL;
  mode = BIG;
  for (i=0; i<ncat; i++)
    {
    icat = catindex + i;
    if (ext == ALL_EXTENSIONS)
      for (e=0; e<next; e++)
        set = read_samples(set, filename[icat], fwhmmin[i]/2.0, fwhmmax[i]/2.0,
			e, next, icat, context, context->pc+i*context->npc);
    else
      set = read_samples(set, filename[icat], fwhmmin[i]/2.0, fwhmmax[i]/2.0,
			ext, next, icat, context, context->pc+i*context->npc);
    if (fwhmmode[i]<mode)
      mode = fwhmmode[i];
    }

  set->fwhm = mode;

  sprintf(str, "%d samples loaded.", set->nsample);
//  NFPRINTF(OUTPUT, str);

  if (!set->nsample)
    warning("No appropriate source found!!","");

  free(fwhmmin);
  free(fwhmmax);
  free(fwhmmode);
/*
  if (set->badflags)
    printf("%d detections discarded with bad SExtractor flags\n",
	set->badflags);
  if (set->badsn)
    printf("%d detections discarded with S/N too low\n",
	set->badsn);
  if (set->badfrmin)
    printf("%d detections discarded with FWHM too small\n",
	set->badfrmin);
  if (set->badfrmax)
    printf("%d detections discarded with FWHM too large\n",
	set->badfrmax);
  if (set->badelong)
    printf("%d detections discarded with elongation too high\n",
	set->badelong);
  if (set->badpix)
    printf("%d detections discarded with too many bad pixels\n\n\n",
	set->badpix);
*/
  return set;
  }


/****** compute_fwhmrange *****************************************************
PROTO   float compute_fwhmrange(float *fwhm, int nfwhm,
		float minin, float maxin, float *minout, float *maxout)
PURPOSE Compute the FWHM range associated to a series of FWHM measurements.
INPUT   Pointer to an array of FWHMs,
	number of FWHMs,
	maximum allowed FWHM variation,
	minimum allowed FWHM,
	maximum allowed FWHM,
	pointer to the lower FWHM range (output),
	pointer to the upper FWHM range (output).
OUTPUT  FWHM mode.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 20/03/2008
*/
static float	compute_fwhmrange(float *fwhm, int nfwhm, float maxvar,
		float minin, float maxin, float *minout, float *maxout)
  {
   float	*fwhmt,*fwhmt2,
		df, dfmin,fmin;
   int		i, nw;

/* Sort FWHMs */
   fqmedian(fwhm, nfwhm);

/* Find the mode */
   nw = nfwhm/4;
   if (nw<4)
     nw = 1;
  dfmin = BIG;
  fmin = 0.0;
  fwhmt = fwhm;
  fwhmt2 = fwhm+nw;
  for (i=nfwhm-nw; i--; fwhmt++,fwhmt2++)
    {  
    if ((df = *fwhmt2 - *fwhmt) < dfmin)
      {
      dfmin = df;
      fmin = (*fwhmt2 + *fwhmt)/2.0;
      }
    }

  if (nfwhm<2)
    fmin = *fwhm;

  dfmin = (float)pow((double)maxvar+1.0, 0.3333333);
  *minout = dfmin>0.0?fmin/dfmin:0.0;
  if (*minout<minin)
    *minout = minin;
  *maxout = fmin*dfmin*dfmin;
  if (*maxout>maxin)
    *maxout = maxin;

  return fmin;
  }


/****** read_samples *********************************************************
PROTO	setstruct *read_samples(setstruct *set, char *filename,
			float frmin, float frmax,
			int ext, int next, int catindex,
			contextstruct *context, double *pcval)

PURPOSE	Read point source data for a given set.
INPUT	Pointer to the data set,
	catalogue filename,
	minimum flux radius,
	maximum flux radius,
	current extension,
	number of extensions,
	catalogue index
	pointer to the context,
	pointer to the array of principal components, if available.
OUTPUT  Pointer to a set containing samples that match acceptance criteria.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/02/2014
*/
setstruct *read_samples(setstruct *set, char *filename,
			float frmin, float frmax,
			int ext, int next, int catindex,
			contextstruct *context, double *pcval)

  {
   catstruct		*cat;
   tabstruct		*tab, *keytab;
   keystruct		*key, *vigkey;
   samplestruct		*sample;
   t_type		contexttyp[MAXCONTEXT];
   void			*contextvalp[MAXCONTEXT];
   static char		str[MAXCHAR], str2[MAXCHAR];
   char			rkeyname[MAXCHAR],
			**kstr,
			*head, *buf, *pstr;
   unsigned int		*imaflags;
   unsigned short	*flags, *wflags;
   double		contextval[MAXCONTEXT],
			*dxm,*dym, *cmin, *cmax, dval;
   float		*xm, *ym, *vignet,*vignett, *flux, *fluxrad, *elong,
			*snr,
			backnoise, backnoise2, gain, minsn,maxelong;
   static int		ncat;
   int			*lxm,*lym,
			i,j, n, nsample,nsamplemax,
			vigw, vigh, vigsize, nobj, nt, ndet,
			maxbad, maxbadflag, ldflag, ext2, pc, contflag;
   short 		*sxm,*sym;


  maxbad = prefs.badpix_nmax;
  maxbadflag = prefs.badpix_flag;
  maxelong = (float)(prefs.maxellip < 1.0?
	(prefs.maxellip + 1.0)/(1.0 - prefs.maxellip)
	: 100.0);
  minsn = prefs.minsn;

/* If a NULL pointer is provided, we allocate a new set */
  if (!set)
    {
    set = init_set(context);
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
  dxm = dym = NULL;
  xm = ym = NULL;
  lxm = lym = NULL;
  sxm = sym = NULL;
  keytab = init_readobj(tab, &buf);

  if (!(key = name_to_key(keytab, prefs.center_key[0])))
    {
    sprintf(str, "*Error*: %s parameter not found in catalogue ",
	prefs.center_key[0]);
    error(EXIT_FAILURE, str, filename);
    }
  if (key->ttype == T_DOUBLE)
    dxm = (double *)key->ptr;
  else if (key->ttype == T_FLOAT)
    xm = (float *)key->ptr;
  else if (key->ttype == T_LONG)
    lxm = (int *)key->ptr;
  else
    sxm = (short *)key->ptr;

  if (!(key = name_to_key(keytab, prefs.center_key[1])))
    {
    sprintf(str, "*Error*: %s parameter not found in catalogue ",
	prefs.center_key[0]);
    error(EXIT_FAILURE, str, filename);
    }
   if (key->ttype == T_DOUBLE)
    dym = (double *)key->ptr;
  else if (key->ttype == T_FLOAT)
    ym = (float *)key->ptr;
  else if (key->ttype == T_LONG)
    lym = (int *)key->ptr;
  else
    sym = (short *)key->ptr;

  if (!(key = name_to_key(keytab, "FLUX_RADIUS")))
    error(EXIT_FAILURE, "*Error*: FLUX_RADIUS parameter not found in catalog ",
		filename);
  fluxrad = (float *)key->ptr;

  strcpy(rkeyname, prefs.photflux_key);
  strtok(rkeyname, "([{}])");
  n = (pstr = strtok(NULL,"([{}])"))? atoi(pstr) - 1 : 0;
  if (!(key = name_to_key(keytab, rkeyname)))
    {
    sprintf(str, "*Error*: %s parameter not found in catalogue ",
	prefs.photflux_rkey);
    error(EXIT_FAILURE, str, filename);
    }
  flux = (float *)key->ptr;
  if (n)
    {
    if (key->naxis==1 && n<key->naxisn[0])
      flux += n;
    else
      {
      sprintf(str, "Not enough apertures for %s in catalogue %s: ",
	prefs.photflux_key, filename);
      warning(str, "using first aperture instead");
      }
    }

  if (!(key = name_to_key(keytab, "SNR_WIN")))
    {
    sprintf(str, "*Error*: SNR_WIN parameter not found in catalogue ");
    error(EXIT_FAILURE, str, filename);
    }
  snr = (float *)key->ptr;

  if ((key = name_to_key(keytab, "ELONGATION")))
    elong = (float *)key->ptr;
  else
    elong = NULL;

/* Load optional SExtractor FLAGS parameter */
  if ((key = name_to_key(keytab, "FLAGS")))
    flags = (unsigned short *)key->ptr;
  else
    flags = NULL;

/* Load optional SExtractor FLAGS_WEIGHT parameter */
  if ((key = name_to_key(keytab, "FLAGS_WEIGHT")))
    wflags = (unsigned short *)key->ptr;
  else
    wflags = NULL;

/* Load optional SExtractor IMAFLAGS_ISO parameter */
  if ((key = name_to_key(keytab, "IMAFLAGS_ISO")))
    imaflags = (unsigned int *)key->ptr;
  else
    imaflags = NULL;

  if (!(key = name_to_key(keytab, "VIGNET")))
    error(EXIT_FAILURE,
	"*Error*: VIGNET parameter not found in catalog ", filename);
  vignet = (float *)key->ptr;
  nobj = key->nobj;

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
  kstr = context->name;
  pc = 0;
  for (i=0; i<set->ncontext; i++, kstr++)
    if (context->pcflag[i])
      {
      contextvalp[i] = &pcval[pc++];
      contexttyp[i] = T_DOUBLE;
      }
    else if (**kstr==(char)':')
      {
      contextvalp[i] = &contextval[i];
      contexttyp[i] = T_DOUBLE;
      if (fitsread(head, *kstr+1, contextvalp[i], H_FLOAT,T_DOUBLE)
		== RETURN_ERROR)
        {
        sprintf(str, "*Error*: %s parameter not found in the header of ",
		*kstr+1);
        error(EXIT_FAILURE, str, filename);
        }
      }
    else
      {
      strcpy(rkeyname, *kstr);
      strtok(rkeyname,"([{}])");
      n= (pstr = strtok(NULL,"([{}])"))? atoi(pstr) - 1 : -1;
      if (!(key = name_to_key(keytab, rkeyname)))
        {
        sprintf(str, "*Error*: %s parameter not found in catalog ", *kstr);
        error(EXIT_FAILURE, str, filename);
        }
      contextvalp[i] = key->ptr;
      contexttyp[i] = key->ttype;
      if (n >= 0)
        {
        if (key->naxis==1 && n<key->naxisn[0])
          contextvalp[i] += n * (key->nbytes / key->naxisn[0]);
        else
          {
          sprintf(str, "Not enough %s elements in catalogue %s: ", rkeyname,
		filename);
          warning(str, "using first element instead");
          }
        }
      strcpy(set->contextname[i], key->name);
      }
  if (next>1)
    sprintf(str2, "[%d/%d]", ext+1, next);
  else
    strcpy(str2, "");

/* Now examine each vector of the shipment */
  nt = keytab->naxisn[1];
  ndet = 0;
  for (n=0; nt; n++)
    {
    ndet++;
    nt = read_obj(keytab,tab, buf);
    if (!(n%100))
      {
      sprintf(str,"Catalog #%d %s: Object #%d / %d samples stored",
	ncat, str2, n,nsample);
//      NFPRINTF(OUTPUT, str);
      }

/*---- Apply some selection over flags, fluxes... */
    contflag = 0;
    if (flags && (*flags&prefs.flag_mask))
      {
      contflag++;
      set->badflags++;
      }
    if (wflags && (*wflags&prefs.wflag_mask))
      {
      contflag++;
      set->badwflags++;
      }
    if (imaflags && (*imaflags&prefs.imaflag_mask))
      {
      contflag++;
      set->badwflags++;
      }
    if (*snr<minsn)
      {
      contflag++;
      set->badsn++;
      }
    if (*fluxrad<frmin)
      {
      contflag++;
      set->badfrmin++;
      }
    if (*fluxrad>frmax)
      {
      contflag++;
      set->badfrmax++;
      }
    if (elong && *elong>maxelong)
      {
      contflag++;
      set->badelong++;
      }
    if (contflag)
      continue;
/*-- ... and check the integrity of the sample */
    j = 0;
    vignett = vignet;
    for (i=vigsize; i--; vignett++)
      if (*vignett <= -BIG)
        j++;
    if (maxbadflag && j > maxbad)
      {
      set->badpix++;
      continue; 
      }
    
/*-- Allocate memory for the first shipment */
    if (!set->sample)
      {
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

/*-- Increase storage space to receive new candidates if needed */
    if (nsample>=nsamplemax)
      {
       int	nadd=(int)(1.62*nsamplemax);
      nsamplemax = nadd>nsamplemax?nadd:nsamplemax+1;
      realloc_samples(set, nsamplemax);
      }

    sample = set->sample + nsample;
    sample->detindex = ndet;
    sample->extindex = ext;
    sample->catindex = catindex;

/*-- Copy the vignet to the training set */
    memcpy(sample->vig, vignet, vigsize*sizeof(float));

    sample->norm = *flux;
    sample->backnoise2 = backnoise2;
    sample->gain = gain;
    if (dxm)
      sample->x = *dxm;
    else if (xm)
      sample->x = *xm;
    else if (lxm)
      sample->x = *lxm;
    else
      sample->x = *sxm;
    if (dym)
      sample->y = *dym;
    else if (ym)
      sample->y = *ym;
    else if (lym)
      sample->y = *lym;
    else
      sample->y = *sym;
    sample->dx = sample->x - (int)(sample->x+0.49999);
    sample->dy = sample->y - (int)(sample->y+0.49999);
    for (i=0; i<set->ncontext; i++)
      {
      dval = sample->context[i];
      ttypeconv(contextvalp[i], &dval, contexttyp[i], T_DOUBLE);
      sample->context[i] = dval;
/*---- Update min and max */
      if (dval<cmin[i])
        cmin[i] = dval;
      if (dval>cmax[i])
        cmax[i] = dval;
      }
    make_weights(set, sample);
    recenter_sample(sample, set, *fluxrad);
    nsample++;
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


/****** recenter_sample ******************************************************
PROTO   void recenter_samples(samplestruct sample,
		setstruct *set, float fluxrad)
PURPOSE Recenter sample image using windowed barycenter.
INPUT   sample structure pointer,
	set structure pointer,
	flux radius.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 16/11/2009
*/
void	recenter_sample(samplestruct *sample, setstruct *set, float fluxrad)

  {
   double	tv, dxpos, dypos;
   float	*ima,*imat,*weight,*weightt,
		pix, var, locpix, locarea, sig,twosig2, raper,raper2,
		offsetx,offsety, mx,my, mx2ph,my2ph,
		rintlim,rintlim2,rextlim2, scalex,scaley,scale2,
		dx,dy, dx1,dy2, r2;
   int		i, x,y,x2,y2, w,h, sx,sy, xmin,xmax,ymin,ymax, pos;

  sig = fluxrad*2.0/2.35; /* From half-FWHM to sigma */
  twosig2 = 2.0*sig*sig;

  w = set->vigsize[0];
  h = set->vigsize[1];
/* Integration radius */
  raper = RECENTER_NSIG*sig;
  raper2 = raper*raper;
/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
  rintlim = raper - 0.75;
  rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
  rextlim2 = (raper + 0.75)*(raper + 0.75);
  scaley = scalex = 1.0/RECENTER_OVERSAMP;
  scale2 = scalex*scaley;
  offsetx = 0.5*(scalex-1.0);
  offsety = 0.5*(scaley-1.0);
/* Use isophotal centroid as a first guess */
  mx = sample->dx + (float)(w/2);
  my = sample->dy + (float)(h/2);

  for (i=0; i<RECENTER_NITERMAX; i++)
    {
    xmin = (int)(mx-raper+0.499999);
    xmax = (int)(mx+raper+1.499999);
    ymin = (int)(my-raper+0.499999);
    ymax = (int)(my+raper+1.499999);
    mx2ph = mx*2.0 + 0.49999;
    my2ph = my*2.0 + 0.49999;

    if (xmin < 0)
      xmin = 0;
    if (xmax > w)
      xmax = w;
    if (ymin < 0)
      ymin = 0;
    if (ymax > h)
      ymax = h;
    tv = 0.0;
    dxpos = dypos = 0.0;
    ima = sample->vig;
    weight = sample->vigweight;
    for (y=ymin; y<ymax; y++)
      {
      imat= ima + (pos = y*w + xmin);
      weightt = weight + pos;
      dy = y - my;
      for (x=xmin; x<xmax; x++, imat++, weightt++)
        {
        dx = x - mx;
        if ((r2=dx*dx+dy*dy)<rextlim2)
          {
          if (RECENTER_OVERSAMP>1 && r2> rintlim2)
            {
            dx += offsetx;
            dy += offsety;
            locarea = 0.0;
            for (sy=RECENTER_OVERSAMP; sy--; dy+=scaley)
              {
              dx1 = dx;
              dy2 = dy*dy;
              for (sx=RECENTER_OVERSAMP; sx--; dx1+=scalex)
                if (dx1*dx1+dy2<raper2)
                  locarea += scale2;
              }
            }
          else
            locarea = 1.0;
          locarea *= expf(-r2/twosig2);
/*-------- Here begin tests for pixel and/or weight overflows. Things are a */
/*-------- bit intricated to have it running as fast as possible in the most */
/*-------- common cases */
          pix = *imat;
          if ((var=*weightt)==0.0)
            {
            if ((x2=(int)(mx2ph-x))>=0 && x2<w
		&& (y2=(int)(my2ph-y))>=0 && y2<h
		&& (var=*(weight + (pos = y2*w + x2)))>0.0)
              {
              var = 1.0/var;
              pix = *(ima+pos);
              }
            else
              pix = var = 0.0;
            }
          dx = x - mx;
          dy = y - my;
          locpix = locarea*pix;
          tv += locpix;
          dxpos += locpix*dx;
          dypos += locpix*dy;
          }
        }
      }

    if (tv>0.0)
      {
      mx += (dxpos /= tv)*RECENTER_GRADFAC;
      my += (dypos /= tv)*RECENTER_GRADFAC;
      }
    else
      break;

/*-- Stop here if position does not change */
    if (dxpos*dxpos+dypos*dypos < RECENTER_STEPMIN*RECENTER_STEPMIN)
      break;
    }

  sample->dx = mx - (float)(w/2);
  sample->dy = my - (float)(h/2);

  return;
  }


/****** malloc_samples *******************************************************
PROTO   void malloc_samples(setstruct *set, int nsample)
PURPOSE Allocate memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
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
    QMALLOC(sample->vigchi, float, set->nvig);
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
VERSION 26/03/2008
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
      QMALLOC(sample->vigchi, float, set->nvig);
      QMALLOC(sample->vigweight, float, set->nvig);
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
      free(sample->vigchi);
      free(sample->vigweight);
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
VERSION 26/03/2008
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
    free(sample->vigchi);
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
VERSION 26/03/2008
*/
setstruct	*init_set(contextstruct *context)

  {
   setstruct	*set;
   int		i;

  QCALLOC(set, setstruct, 1);
  set->nsample = set->nsamplemax = 0;
  set->vigdim = 2;
  QMALLOC(set->vigsize, int, set->vigdim);
  set->ncontext = context->ncontext;
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
VERSION 26/03/2008
*/
void	end_set(setstruct *set)

  {
   int	i;

  free_samples(set);
  free(set->vigsize);
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

