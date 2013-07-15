/*
*				field.c
*
* Manage multiple PSFs.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		25/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
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
#include	"fitswcs.h"
#include	"misc.h"
#include	"prefs.h"
#include	"psf.h"
#include	"field.h"


/****** field_init_finalize ************************************************************
PROTO	void field_init(fieldstruct *)
PURPOSE	Finish allocating and initializing a PSF MEF structure (groups of PSFs).
INPUT	fieldstruct pointer
OUTPUT  -.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 08/04/2010
 ***/
void
field_init_finalize(fieldstruct *field)
{
  int countsize = prefs.context_nsnap*prefs.context_nsnap;
  int e, next0 = field->next;

  field_locate(field);
  QCALLOC(field->ccat, catstruct *, MAXCHECK);
  countsize = prefs.context_nsnap*prefs.context_nsnap;
  QMALLOC(field->lcount, int *, next0);
  QMALLOC(field->acount, int *, next0);
  QMALLOC(field->count, int *, next0);
  QMALLOC(field->modchi2, double *, next0);
  QMALLOC(field->modresi, double *, next0);
  for (e=0; e<next0; e++)
    {
    QCALLOC(field->lcount[e], int, countsize);
    QCALLOC(field->acount[e], int, countsize);
    QCALLOC(field->count[e], int, countsize);
    QCALLOC(field->modchi2[e], double, countsize);
    QCALLOC(field->modresi[e], double, countsize);
    }
  }
