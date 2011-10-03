/*
*				check.h
*
* Include file for check.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		03/10/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef	_FIELD_H_
#include "field.h"
#endif

#ifndef _CHECK_H_
#define _CHECK_H_

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHECK	16		/* max. # of CHECKimages */

/*----------------------------- Type definitions --------------------------*/
typedef enum {PSF_NONE, PSF_BASIS, PSF_CHI, PSF_PROTO, PSF_RESIDUALS,
		PSF_RESIDUALS_GRID, PSF_SAMPLES, PSF_SAMPLES_GRID,
		PSF_SNAPSHOTS, PSF_SNAPSHOTS_IMRES,
		PSF_WEIGHTS, PSF_MOFFAT,PSF_SUBMOFFAT,PSF_SUBSYM, PSF_GREAT}
	checkenum;

/*---------------------------------- protos --------------------------------*/
extern void		check_write(fieldstruct *field,	setstruct *set,
				char *checkname, checkenum checktype,
				int ext, int next, int cubeflag);

#endif

