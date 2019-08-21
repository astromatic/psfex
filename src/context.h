/*
*				context.h
*
* Include file for context.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2019 IAP/CNRS/SorbonneU
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
*	Last modified:		21/08/2019
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*--------------------------------- constants -------------------------------*/

#define		MAXCONTEXT		8	/* max. # of context keys */
#define		CONTEXT_KEEPHIDDEN	0
#define		CONTEXT_REMOVEHIDDEN	1

/*--------------------------- structure definitions -------------------------*/

typedef struct context
  {
  char		**name;			/* Context names */
  int		*group;			/* Context groups */
  int		*pcflag;		/* Flags PC contexts */
  int		ncontext;		/* Total number of contexts */
  int		*degree;		/* Group degrees */
  int		ngroup;			/* Number of context groups */
  double	*pc;			/* PC components */
  int		npc;			/* Number of PC components */
  }	contextstruct;

typedef struct psf psfstruct;

/*-------------------------------- protos -----------------------------------*/

contextstruct	*context_init(char **names, int *group, int ndim, int *degree,
			int ngroup, int pcexflag);
void 		context_end(contextstruct *context),
		context_to_pos(double *context, wcsstruct *wcs,
			psfstruct *psf, double *pos),
		pos_to_context(double *pos, psfstruct *psf,
			wcsstruct *wcs, double *context);

#endif

