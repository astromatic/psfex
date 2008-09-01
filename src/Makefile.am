# Program Makefile for PSFEx
# Copyright (C) 2002-2008 Emmanuel Bertin.
SUBDIRS			= fits levmar wcs
if USE_PLPLOT
CPLOTSOURCE		= cplot.c
endif
bin_PROGRAMS		= psfex
psfex_SOURCES		= check.c context.c $(CPLOTSOURCE) diagnostic.c fft.c \
			  field.c fitswcs.c homo.c main.c makeit.c misc.c \
			  pca.c poly.c prefs.c psf.c sample.c vignet.c xml.c \
			  check.h context.h cplot.h define.h diagnostic.h \
			  fft.h field.h fitswcs.h globals.h homo.h key.h pca.h \
			  poly.h prefs.h preflist.h psf.h sample.h threads.h \
			  types.h vignet.h wcscelsys.h xml.h
psfex_LDADD		= $(top_builddir)/src/fits/libfits.a \
			  $(top_builddir)/src/levmar/liblevmar.a \
			  $(top_builddir)/src/wcs/libwcs_c.a

DATE=`date +"%Y-%m-%d"`
