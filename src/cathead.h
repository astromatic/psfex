/*
*				cathead.h
*
* Merged and full catalogue headers.
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

#ifndef _CATOUT_H_
#include "catout.h"
#endif

/* Output catalog fields */

outsamplestruct	refoutsample;
keystruct	refoutkey[] = {
  {"SOURCE_NUMBER", "Source index",
	&refoutsample.detindex, H_INT, T_LONG,
	"%10d", "", "meta.number", ""},
   {"EXTENSION", "Extension index",
	&refoutsample.extindex, H_INT, T_SHORT,
	"%4d", "", "meta.number", ""},
   {"CATALOG_NUMBER", "File index",
	&refoutsample.catindex, H_INT, T_LONG,
	"%7d", "", "meta.number", ""},
  {"VECTOR_CONTEXT", "Context vector",
	&refoutsample.context, H_FLOAT, T_DOUBLE,
	"%12.6g", "", "obs.param", "",
	1, &refoutsample.ncontext},
   {"X_IMAGE", "Position along x image axis",
	&refoutsample.x, H_FLOAT, T_DOUBLE,
	"%11.4f", "pixel", "pos.cartesian.x", "pix"},
   {"Y_IMAGE", "Position along y image axis",
	&refoutsample.y, H_FLOAT, T_DOUBLE,
	"%11.4f", "pixel", "pos.cartesian.y", "pix"},
   {"DELTAX_IMAGE", "Position offset along x image axis",
	&refoutsample.dx, H_FLOAT, T_FLOAT,
	"%11.4f", "pixel", "pos.cartesian.x;arith.diff", "pix"},
   {"DELTAY_IMAGE", "Position offset along y image axis",
	&refoutsample.dy, H_FLOAT, T_FLOAT,
	"%11.4f", "pixel", "pos.cartesian.y;arith.diff", "pix"},
   {"NORM_PSF", "Source (inverse) normalization factor",
	&refoutsample.norm, H_FLOAT, T_FLOAT,
	"%12.6g", "count", "phot.flux;instr.det.psf", "ct"},
   {"CHI2_PSF", "PSF fitting chi2/d.o.f.",
	&refoutsample.chi2, H_FLOAT, T_FLOAT,
	"%12.6g", "", "stat.fit.chi2;instr.det.psf", ""},
   {"RESI_PSF", "PSF fitting normalized residuals",
	&refoutsample.modresi, H_FLOAT, T_FLOAT,
	"%12.6g", "", "stat.fit.residual;instr.det.psf", ""},
  {""},
  };

