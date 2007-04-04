# Program Makefile for PSFEx
# Copyright (C) 2002-2007 Emmanuel Bertin.
SUBDIRS			= fits lmfit levmar
bin_PROGRAMS		= psfex
psfex_SOURCES		= check.c diagnostic.c fft.c main.c makeit.c misc.c \
			  pc.c poly.c prefs.c psf.c sample.c vignet.c xml.c \
			  check.h define.h diagnostic.h fft.h globals.h key.h \
			  poly.h prefs.h preflist.h psf.h sample.h types.h \
			  vignet.h xml.h
AM_CPPFLAGS		= -DXSL_URL=\"file://$(pkgdatadir)/${PACKAGE_NAME}.xsl\"
psfex_LDADD		= $(top_builddir)/src/fits/libfits.a \
			  $(top_builddir)/src/lmfit/liblmfit.a \
			  $(top_builddir)/src/levmar/liblevmar.a
DATE=`date +"%Y-%m-%d"`