/*
*				wcs_utils.c
*
* Utilities to help manage World Coordinate System data.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		20/11/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <math.h>

#include "define.h"
struct tab; typedef struct tab tabstruct;

#include	"fitswcs.h"

/******* wcs_dist ***********************************************************
PROTO	double wcs_dist(int naxis, int lat, int lng, double *wcspos1, double *wcspos2)
PURPOSE	Compute the angular distance between 2 points on the sky.
INPUT	naxis, lat, lng unpacked from a WCS structure,
	Pointer to the first array of world coordinates,
	Pointer to the second array of world coordinates.
OUTPUT	Angular distance (in degrees) between points.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/07/2002
 ***/
double	wcs_dist_impl(int naxis, int lat, int lng, double *wcspos1, double *wcspos2)

  {
  double	d, dp;
  int		i;

  if (lat!=lng)
    {
/*-- We are operating in angular coordinates */
    d = sin(wcspos1[lat]*DEG)*sin(wcspos2[lat]*DEG)
	+ cos(wcspos1[lat]*DEG)*cos(wcspos2[lat]*DEG)
		*cos((wcspos1[lng]-wcspos2[lng])*DEG);
    return d>-1.0? (d<1.0 ? acos(d)/DEG : 0.0) : 180.0;
    }
  else
    {
    d = 0.0;
    for (i=0; i<naxis; i++)
      {
      dp = wcspos1[i] - wcspos2[i];
      d += dp*dp;
      }
    return sqrt(d);
    }
  }
