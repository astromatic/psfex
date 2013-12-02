/*
*				globals.h
*
* Global declarations.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		02/12/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CONTEXT_H_
#include "context.h"
#endif
#ifndef _PSF_H_
#include "psf.h"
#endif
#ifndef _SAMPLE_H_
#include "sample.h"
#endif

/*----------------------- miscellaneous variables ---------------------------*/
char		gstr[MAXCHAR];

/*------------------------------- functions ---------------------------------*/
extern  void	makeit(void);
psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context);

