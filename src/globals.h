/*
*				globals.h
*
* Global declarations.
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
*	Last modified:		06/01/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#include	<time.h>

/*----------------------- miscellaneous variables ---------------------------*/
extern char		gstr[MAXCHAR];

/*------------------------------- functions ---------------------------------*/
#include "fits/fitscat.h"
#include "psf.h"
#include "sample.h"
#include "field.h"

extern  void	error(int, const char *, const char *),
		makeit(void),
		makeit_body(fieldstruct **fields, contextstruct **context,
                            contextstruct **fullcontext, int free_sets);
psfstruct	*make_psf(setstruct *set, float psfstep,
                          float *basis, int nbasis, contextstruct *context);
