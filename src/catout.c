/*
*				catout.c
*
* Write catalogs containing the set of stars used for computing the PSF model.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	along with PSFEx. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		26/02/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "cathead.h"
#include "catout.h"
#include "prefs.h"
#include "context.h"
#include "sample.h"

/****** init_outcat *************************************************
PROTO	outcatstruct *init_outcat(char *filename, int ncontext)
PURPOSE	Initialize a SExtractor-like catalog containing detections used by PSFex
	for computing the PSF model.
INPUT	File name,
	number of contexts.
OUTPUT  Pointer to the output catalog, ready for writing.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 26/02/2014
*/
outcatstruct	*init_outcat(char *filename, int ncontext)

  {
   static char  imtabtemplate[][80] = {
"SIMPLE  =                    T / This is a FITS file",
"BITPIX  =                    8 / ",
"NAXIS   =                    2 / 2D data",
"NAXIS1  =                    1 / Number of rows",
"NAXIS2  =                    1 / Number of columns",
"EXTEND  =                    T / This file may contain FITS extensions",
"END                            "};
   catstruct		*cat;
   tabstruct		*asctab, *imtab, *objtab;
   keystruct		*key, *objkeys;
   samplestruct		*samp;
   outsamplestruct	*outsample;
   outcatstruct		*outcat;
   FILE			*ascfile;
   char			str[80],
			*buf, *rfilename;
   long			dptr;
   int			i,k,n;

  if (prefs.outcat_type == CAT_NONE)
    return NULL;

  QCALLOC(outcat, outcatstruct, sizeof(outcatstruct));
  outsample = &outcat->outsample;
  outcat->ncontext = refoutsample.ncontext = ncontext;

/* LDAC Object header */
  outcat->objtab = objtab = new_tab("LDAC_OBJECTS");
/* Set key pointers */
  QCALLOC(objkeys, keystruct, (sizeof(refoutkey) / sizeof(keystruct)));
  dptr = (long)((char *)outsample - (char *)&refoutsample);
  for (k=0; refoutkey[k].name[0]; k++)
    {
    objkeys[k] = refoutkey[k];
    key = objkeys+k;
/*-- A trick to access the fields of the dynamic outsample structure */
    key->ptr = (void *)((char *)key->ptr + dptr);
    key->nbytes = t_size[key->ttype]*(key->naxis? *key->naxisn : 1);
    add_key(key,objtab, 0);
    }
/* Create a new output catalog */
  if (prefs.outcat_type == CAT_ASCII_HEAD
	|| prefs.outcat_type == CAT_ASCII
	|| prefs.outcat_type == CAT_ASCII_VOTABLE)
    {
    cat = NULL;
    if (prefs.outcatpipe_flag)
      ascfile = stdout;
    else
      if (!(ascfile = fopen(filename, "w+")))
        error(EXIT_FAILURE,"*Error*: cannot open ", filename);
    if (prefs.outcat_type == CAT_ASCII_HEAD && (key = objtab->key))
      for (i=0,n=1; i++<objtab->nkey; key=key->nextkey)
        {
        if (*key->unit)
          fprintf(ascfile, "# %3d %-22.22s %-58.58s [%s]\n",
                n, key->name,key->comment, key->unit);
        else
          fprintf(ascfile, "# %3d %-22.22s %-58.58s\n",
                n, key->name,key->comment);
        n += key->naxis? *key->naxisn : 1;
        }
    else if (prefs.outcat_type == CAT_ASCII_VOTABLE && objtab->key) 
      {
/*---- A short, "relative" version of the filename */
      if (!(rfilename = strrchr(filename, '/')))
        rfilename = filename;
      else
        rfilename++;
      write_xml_header(ascfile);
      fprintf(ascfile,
	" <TABLE ID=\"Output_List\" name=\"%s/out\">\n", rfilename);
      fprintf(ascfile,
        "  <DESCRIPTION>Table of detections used by %s</DESCRIPTION>\n",
	BANNER);
      fprintf(ascfile,
        "  <!-- Now comes the definition of each %s parameter -->\n", BANNER);
      write_vo_fields(ascfile, objtab);
      fprintf(ascfile, "   <DATA><TABLEDATA>\n");
      }
    }
  else
    {
    ascfile = NULL;
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);

/*-- Primary header */
    save_tab(cat, cat->tab);

/*-- We create a dummy table (only used through its header) */
    QCALLOC(asctab, tabstruct, 1);
    asctab->headnblock = 1 + (sizeof(imtabtemplate)-1)/FBSIZE;
    QCALLOC(asctab->headbuf, char, asctab->headnblock*FBSIZE);
    memcpy(asctab->headbuf, imtabtemplate, sizeof(imtabtemplate));
    for (buf = asctab->headbuf, i=FBSIZE*asctab->headnblock; i--; buf++)
      if (!*buf)
        *buf = ' ';
/*-- (dummy) LDAC Image header */

    imtab = new_tab("LDAC_IMHEAD");
    key = new_key("Field Header Card");
    key->ptr = asctab->headbuf;
    asctab->headbuf = NULL;
    free_tab(asctab);
    key->naxis = 2;
    QMALLOC(key->naxisn, int, key->naxis);
    key->naxisn[0] = 80;
    key->naxisn[1] = fitsfind(key->ptr, "END     ")+1;
    key->htype = H_STRING;
    key->ttype = T_STRING;
    key->nobj = 1;
    key->nbytes = key->naxisn[0]*key->naxisn[1];
    add_key(key, imtab, 0);
    save_tab(cat, imtab);
    free_tab(imtab);
    objtab->cat = cat;
    init_writeobj(cat, objtab, &buf);
    }

  outcat->ascfile = ascfile;
  outcat->objkeys = objkeys;
  outcat->buf = buf;

  return outcat;
  }


/****** write_outcat *************************************************
PROTO	void write_outcat(outcatstruct *outcat, setstruct *set)
PURPOSE	Write the content of a set in the output catalog.
INPUT	Pointer to the output catalog structure,
	pointer to the catalogue set.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 26/02/2014
*/
void	write_outcat(outcatstruct *outcat, setstruct *set)

  {
   outsamplestruct	*outsample;
   samplestruct		*samp;
   int			c, n, nc;

  nc = outcat->ncontext;
  outsample = &outcat->outsample;
  samp = set->sample;
  for (n=set->nsample; n--; samp++)
    {
    memset(outsample, 0, sizeof(outsample));
    outsample->detindex = samp->detindex;
    outsample->extindex = samp->extindex + 1;
    outsample->catindex = samp->catindex + 1;
    for (c=0; c<nc; c++)
      outsample->context[c] = samp->context[c];
    outsample->x = samp->x;
    outsample->y = samp->y;
    outsample->dx = samp->dx;
    outsample->dy = samp->dy;
    outsample->norm = samp->norm;
    outsample->chi2 = samp->chi2;
    outsample->modresi = samp->modresi; 

/*-- Write to the catalog */
    if (prefs.outcat_type == CAT_ASCII_HEAD || prefs.outcat_type == CAT_ASCII)
      print_obj(outcat->ascfile, outcat->objtab);
    else if (prefs.outcat_type == CAT_ASCII_VOTABLE)
      voprint_obj(outcat->ascfile, outcat->objtab);
    else
      write_obj(outcat->objtab, outcat->buf);
    }

  return;
  }


/****** end_outcat *************************************************
PROTO	void end_outcat(outcatstruct *outcat)
PURPOSE	End output catalog.
INPUT	Pointer to the output catalog structure.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 26/02/2014
*/
void	end_outcat(outcatstruct *outcat)

  {
   tabstruct	*objtab;

  objtab = outcat->objtab;
  if (prefs.outcat_type == CAT_ASCII_HEAD || prefs.outcat_type == CAT_ASCII)
    {
    if (!prefs.outcatpipe_flag)
      fclose(outcat->ascfile);
    }
  else if (prefs.outcat_type == CAT_ASCII_VOTABLE)
    {
    fprintf(outcat->ascfile, "    </TABLEDATA></DATA>\n");
    fprintf(outcat->ascfile, "  </TABLE>\n");
/*-- Add configuration file meta-data */
    write_xml_meta(outcat->ascfile, NULL);
    fprintf(outcat->ascfile, "</RESOURCE>\n");
    fprintf(outcat->ascfile, "</VOTABLE>\n");
    }
  else
    end_writeobj(objtab->cat, objtab, outcat->buf);

  objtab->key = NULL;
  objtab->nkey = 0;
  free(outcat->objkeys);

  if (objtab->cat)
    free_cat(&objtab->cat, 1);
  else
    free_tab(objtab);

  free(outcat);

  return;
  }


/****** write_vo_fields *******************************************************
PROTO	void	write_vo_fields(FILE *file, tabstruct *objtab)
PURPOSE	Write the list of columns to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream),
	Pointer to the object table.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	19/10/2009
 ***/
void	write_vo_fields(FILE *file, tabstruct *objtab)
  {
   keystruct	*key;
   char		datatype[40], arraysize[40], str[40];
   int		i, d;

  if (!objtab || !objtab->key)
    return;
  key=objtab->key;
  for (i=0; i++<objtab->nkey; key=key->nextkey)
    {
/*--- indicate datatype, arraysize, width and precision attributes */
/*--- Handle multidimensional arrays */
    arraysize[0] = '\0';
    if (key->naxis>1)
      {
      for (d=0; d<key->naxis; d++)
        {
        sprintf(str, "%s%d", d?"x":" arraysize=\"", key->naxisn[d]);
        strcat(arraysize, str);
        }
      strcat(arraysize, "\"");
      }
    switch(key->ttype)
      {
      case T_BYTE:	strcpy(datatype, "unsignedByte"); break;
      case T_SHORT:	strcpy(datatype, "short"); break;
      case T_LONG:	strcpy(datatype, "int"); break;
      case T_FLOAT:	strcpy(datatype, "float"); break;
      case T_DOUBLE:	strcpy(datatype, "double"); break;
      default:		error(EXIT_FAILURE,
			"*Internal Error*: Unknown datatype in ",
			"initcat()");
      }
    fprintf(file,
	"  <FIELD name=\"%s\" ucd=\"%s\" datatype=\"%s\" unit=\"%s\"%s>\n",
	key->name, key->voucd, datatype,key->vounit, arraysize);
    fprintf(file, "   <DESCRIPTION>%s</DESCRIPTION>\n", key->comment);
    fprintf(file, "  </FIELD>\n");
    }

  return;
  }


