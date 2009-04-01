# Main Makefile for PSFEx
# Copyright (c) 2002-2009 Emmanuel Bertin IAP-CNRS/Univ.P.&M.Curie
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
dist_pkgdata_DATA	= xsl/psfex.xsl
EXTRA_DIST		= config doc AUTHORS BUGS ChangeLog \
			  COPYRIGHT HISTORY INSTALL README THANKS \
			  acx_atlas.m4 acx_fftw.m4 acx_prog_cc_optim.m4 \
			  acx_plplot.m4 acx_urbi_resolve_dir.m4
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
dist-hook:
	rm -rf `find $(distdir) -name .svn`

rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

