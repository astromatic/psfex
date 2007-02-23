# Program Makefile for PSFEx
# Copyright (C) 2002-2007 Emmanuel Bertin.
SUBDIRS			= fits
bin_PROGRAMS		= psfex
psfex_SOURCES		= check.c fft.c main.c makeit.c misc.c pc.c poly.c \
			  prefs.c psf.c sample.c vignet.c xml.c \
			  check.h define.h fft.h globals.h key.h \
			  poly.h prefs.h preflist.h psf.h sample.h types.h \
			  vignet.h xml.h
psfex_LDADD		= $(top_builddir)/src/fits/libfits.a
DATE=`date +"%Y-%m-%d"`
