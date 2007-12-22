# Main Makefile for PSFEx
# Copyright (C) 2002-2007 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
dist_pkgdata_DATA	= xsl/psfex.xsl
EXTRA_DIST		= config doc AUTHORS BUGS ChangeLog \
			  COPYING HISTORY INSTALL README THANKS \
			  acx_atlas.m4 acx_fftw.m4 acx_prog_cc_optim.m4 \
			  acx_urbi_resolve_dir.m4
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
dist-hook:
	rm -rf `find $(distdir) -name .svn`
rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall \
		-march=opteron" TPXLIBS="-lfftw3 -llapack -lcblas -latlas -lm" \
		rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall \
	-m32 -march=athlon -msse -mfpmath=sse -mtune=athlon" \
	TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -fftw3 -llapack -lcblas \
	-latlas -lm" rpmbuild -ba --target=athlon --clean $(PACKAGE_NAME).spec

rpm-pentium:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall \
	-m32 -march=i686 -msse -mfpmath=sse -mno-3dnow" \
	TPXLIBS="-L/usr/lib/ -L/usr/local/lib/atlas_pentium/ -llapack -lcblas \
	-latlas -L/usr/local/lib/ -lm -lfftw3" rpmbuild -ba --target=i686 \
	--clean $(PACKAGE_NAME).spec

rpm-opteron-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -axWP -ip -unroll" TPXLIBS="-lfftw3 \
	-llapack -lcblas -latlas -lm" rpmbuild -ba --target=x86_64 \
	--clean $(PACKAGE_NAME).spec

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -axKWNPB -ip -unroll" \
	TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lfftw3 -llapack -lcblas \
	-latlas -lm" rpmbuild -ba --target=athlon --clean $(PACKAGE_NAME).spec

rpm-pentium-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -axKWNPB -ip -unroll" \
	TPXLIBS="-L/usr/lib/ -L/usr/local/lib/atlas_pentium/ -llapack -lcblas \
	-latlas -L/usr/local/lib/ -lm -lfftw3" rpmbuild -ba --target=i686 \
	--clean $(PACKAGE_NAME).spec

athlon:
	$(MAKE) CFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall \
	-m32 -march=athlon -msse -mfpmath=sse -mtune=athlon" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lfftw3 -llapack -lcblas -latlas \
	-lm"

pentium:
	$(MAKE) CFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall \
	-m32 -march=i686 -msse -mfpmath=sse -mno-3dnow" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lfftw3 -llapack -lcblas -latlas \
	-lm"

icc:
	$(MAKE) CC="icc" CFLAGS="-O3 -axKNPB -ip -unroll"

athlon-icc:
	$(MAKE) CC="icc32" CFLAGS="-O3 -axKWNPB -ip -unroll" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lfftw3 -llapack -lcblas -latlas \
	-lm"

pentium-icc:
	$(MAKE) CC="icc32" CFLAGS="-O3 -axKWNPB -ip -unroll" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lfftw3 -llapack -lcblas -latlas \
	-lm"

