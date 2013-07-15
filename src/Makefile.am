#
#				Makefile.am
#
# src Makefile.am. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	PSFEx
#
#	Copyright:		(C) 2002-2012 IAP/CNRS/UPMC
#
#	Author:			Emmanuel Bertin (IAP)
#
#	License:		GNU General Public License
#
#	PSFEx is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	PSFEx is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with PSFEx. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		03/04/2013
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBDIRS			= fits levmar wcs
if USE_PLPLOT
CPLOTSOURCE		= cplot.c
endif
bin_PROGRAMS		= psfex
psfex_SOURCES		= check.c context.c $(CPLOTSOURCE) diagnostic.c fft.c \
			  field.c field_utils.c fitswcs.c homo.c main.c makeit.c makeit2.c misc.c \
			  pca.c prefs.c psf.c sample.c sample_utils.c vignet.c wcs_utils.c xml.c \
			  check.h context.h cplot.h define.h diagnostic.h \
			  fft.h field.h fitswcs.h globals.h homo.h key.h \
			  misc.h pca.h prefs.h preflist.h psf.h \
			  sample.h threads.h types.h vignet.h wcscelsys.h xml.h
psfex_LDADD		= $(srcdir)/fits/libfits.a \
			  $(srcdir)/levmar/liblevmar.a \
			  $(srcdir)/wcs/libwcs_c.a

DATE=`date +"%Y-%m-%d"`

