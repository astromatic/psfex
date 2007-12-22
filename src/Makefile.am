# Program Makefile for PSFEx
# Copyright (C) 2002-2007 Emmanuel Bertin.
SUBDIRS			= fits levmar
bin_PROGRAMS		= psfex
psfex_SOURCES		= check.c diagnostic.c fft.c homo.c main.c makeit.c \
			  misc.c pca.c poly.c prefs.c psf.c sample.c vignet.c \
			  xml.c \
			  check.h define.h diagnostic.h fft.h globals.h homo.h \
			  key.h pca.h poly.h prefs.h preflist.h psf.h sample.h \
			  types.h vignet.h xml.h
psfex_LDADD		= $(top_builddir)/src/fits/libfits.a \
			  $(top_builddir)/src/levmar/liblevmar.a
DATE=`date +"%Y-%m-%d"`
