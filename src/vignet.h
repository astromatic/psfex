/*
*				vignet.h
*
* Include file for vignet.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1998-2010 IAP/CNRS/SorbonneU
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
*	Last modified:		22/05/2019
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _VIGNET_H_
#define _VIGNET_H_

/*----------------------------- Internal constants --------------------------*/

#define APER_OVERSAMP	5	/* oversampling in each dimension (MAG_APER) */

#define	INTERP_MAXKERNELWIDTH	8	/* Interpolation kernel buffer size */
#define	INTERPTYPE	INTERP_LANCZOS3	/* Interpolation type */
#define	INTERPW		6	/* Interpolation function range */
#define	INTERPFAC	3.0	/* Interpolation envelope factor */

#define	INTERPF(x)	(x<1e-5 && x>-1e-5? 1.0 \
			:(x>INTERPFAC?0.0:(x<-INTERPFAC?0.0 \
			:sin(PI*x)*sin(PI/INTERPFAC*x)/(PI*PI/INTERPFAC*x*x))))

//#define	INTERPF(x)	(fabs(x)>1.0?0.0 : 1 - fabs(x))
//#define	INTERPF(x)	(fabs(x)>0.5? 0.0:1.0)
/*
#define	INTERPF(x)	(x==0.0?1.0 \
			:(x>INTERPLIM?0.0:(x<-INTERPLIM?0.0 \
			:sin(PI*x)*exp(-x*x/4.4)/(PI*x))))
*/
				/* Lanczos approximation */

/*--------------------------------- typedefs --------------------------------*/

typedef  enum {VIGNET_CPY, VIGNET_ADD, VIGNET_SUB, VIGNET_MUL, VIGNET_DIV}
		vigopenum;

typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;


/*---------------------------------- protos --------------------------------*/
extern void	vignet_make_kernel(float pos, float *kernel,
			interpenum interptype);

extern int	vignet_copy(float *pix1, int *size1,
			float *pix2, int *size2, int idx, int idy,
			vigopenum vigop),
		vignet_resample(float *pix1, int *size1, float *pix2,
			int *size2, double dx, double dy, float step2,
			float stepi, float *dgeoxpix, float *dgeoypix,
			wcsstruct *wcs1, wcsstruct *wcs2);

extern float	vignet_interpolate_pix(float *posin, float *pix, int *naxisn,
			interpenum interptype),
		vignet_aperflux(float *ima, float *var, int w, int h,
			float dxc, float dyc, float aper,
			float gain, float backnoise, float *fluxvar);

#endif

