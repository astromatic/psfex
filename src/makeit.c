/*
*				makeit.c
*
* Main loop.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		19/07/2012
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
#include	"pca.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"xml.h"

time_t		thetime, thetime2;
void		write_error(const char *msg1, const char *msg2);

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

  if (prefs.xml_flag)
    init_xml(ncat);  
  
  makeit_body(fields, &context, &fullcontext, 1);
  next = fields[0]->next;

/* Write XML */
  if (prefs.xml_flag)
    {
    NFPRINTF(OUTPUT, "Writing XML file...");
    write_xml(prefs.xml_name);
    end_xml();
    }

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

/* Free memory */
  for (c=0; c<ncat; c++)
    field_end(fields[c]);
  free(fields);

  if (context->npc)
    context_end(fullcontext);   
  context_end(context);   

  return;
  }




/****** write_error ********************************************************
PROTO	void    write_error(const char *msg1, const char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/02/2007
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
