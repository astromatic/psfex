/*
*				homo.h
*
* Include file for homo.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2008-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _HOMO_H_
#define _HOMO_H_

/*----------------------------- Internal constants --------------------------*/

#define		HOMO_NSNAP	5	/* Number of points per PSFVar dim. */

/*----------------------------- Global variables ---------------------------*/
/*---------------------------------- protos --------------------------------*/
extern void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
			int homobasis_number, double homobasis_scale,
			int ext, int next),
		psf_savehomo(psfstruct *psf, char *filename, int ext, int next);
#endif

