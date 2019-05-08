/*
*				psfbasis.h
*
* Include file for psfbasis.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2019 IAP/CNRS/SorbonneU
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
*	Last modified:		08/05/2019
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _POLY_H_
#include "wcs/poly.h"
#endif

#ifndef _SAMPLE_H_
#include "sample.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _PSFBASIS_H_
#define _PSFBASIS_H_

/*----------------------------- Internal constants --------------------------*/
#define	GAUSS_LAG_OSAMP	3	/* Gauss-Laguerre oversampling factor */


/*----------------------------- Type definitions --------------------------*/
typedef enum {BASIS_NONE, BASIS_PIXEL, BASIS_GAUSS_LAGUERRE, BASIS_FILE,
		BASIS_PIXEL_AUTO}
        basistypenum;
/*--------------------------- structure definitions -------------------------*/


/*---------------------------------- protos --------------------------------*/
extern void	psfbasis_make(psfstruct *psf, setstruct *set,
			basistypenum basis_type,  int nvec);

extern int	psfbasis_pshapelet(float **shape, int w, int h, int nmax,
			double beta),
		psf_readbasis(psfstruct *psf, char *filename, int ext);

#endif

