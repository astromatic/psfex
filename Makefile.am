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

rpm-opteron:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="@PLPLOT_CFLAGS@ -O3 -g -funroll-loops \
		-fomit-frame-pointer -Wall -march=opteron" \
		TPXLDFLAGS="-static -shared-libgcc" \
		TPXLIBS="-lplplotd -lpthread -lfftw3f \
		-llapack -lcblas -latlas -lm" \
		rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="@PLPLOT_CFLAGS@ -O3 -g -funroll-loops \
	-fomit-frame-pointer -Wall -m32 -march=athlon -msse -mfpmath=sse \
	-mtune=athlon" \
	TPXLDFLAGS="-static -shared-libgcc" \
	TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lplplotd -lpthread \
	-fftw3f -llapack -lcblas \
	-latlas -lm" rpmbuild -ba --target=athlon --clean $(PACKAGE_NAME).spec

rpm-pentium:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="@PLPLOT_CFLAGS@ -O3 -g -funroll-loops \
	-fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse \
	-mno-3dnow" \
	TPXLDFLAGS="-static -shared-libgcc" \
	TPXLIBS="-L/usr/lib/ -L/usr/local/lib/atlas_pentium/ \
	-lplplotd -lpthread- llapack -lcblas \
	-latlas -L/usr/local/lib/ -lm -lfftw3f" rpmbuild -ba --target=i686 \
	--clean $(PACKAGE_NAME).spec

rpm-opteron-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="@PLPLOT_CFLAGS@ -O3 -axWPT -ip -unroll" \
	TPXLDFLAGS="-lplplotd -lpthread -static -shared-libgcc -static-intel" \
	TPXLIBS="-lfftw3f \
	-llapack -lcblas -latlas -lm" rpmbuild -ba --target=x86_64 \
	--clean $(PACKAGE_NAME).spec

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="@PLPLOT_CFLAGS@ -O3 -axKWNPT -ip -unroll" \
	TPXLDFLAGS="-static -shared-libgcc -static-intel" \
	TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lplplotd -lpthread \
	-lfftw3f -llapack -lcblas \
	-latlas -lm" rpmbuild -ba --target=athlon --clean $(PACKAGE_NAME).spec

rpm-pentium-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="@PLPLOT_CFLAGS@ -O3 -axKWNPT -ip -unroll" \
	TPXLDFLAGS="-static -shared-libgcc" \
	TPXLIBS="-L/usr/lib/ -L/usr/local/lib/atlas_pentium/ \
	-lplplotd -lpthread -llapack -lcblas \
	-latlas -L/usr/local/lib/ -lm -lfftw3f" rpmbuild -ba --target=i686 \
	--clean $(PACKAGE_NAME).spec

athlon:
	$(MAKE) CFLAGS="@PLPLOT_CFLAGS@ -O3 -g -funroll-loops \
	-fomit-frame-pointer -Wall -m32 -march=athlon -msse -mfpmath=sse \
	-mtune=athlon" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lplplotd -lpthread \
	-lfftw3f -llapack -lcblas -latlas \
	-lm"

pentium:
	$(MAKE) CFLAGS="@PLPLOT_CFLAGS@ -O3 -g -funroll-loops \
	-fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse \
	-mno-3dnow" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lplplotd -lpthread \
	-lfftw3f -llapack -lcblas -latlas \
	-lm"

icc:
	$(MAKE) CC="icc" CFLAGS="@PLPLOT_CFLAGS@ -O3 -axWPT -ip -unroll"

athlon-icc:
	$(MAKE) CC="icc32" CFLAGS="@PLPLOT_CFLAGS@ -O3 -axKWNPT -ip -unroll" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lplplotd -lpthread \
	-lfftw3f -llapack -lcblas -latlas \
	-lm"

pentium-icc:
	$(MAKE) CC="icc32" CFLAGS="@PLPLOT_CFLAGS@ -O3 -axKWNPT -ip -unroll" \
	LIBS="-L/usr/local/lib/ -L/usr/lib/ -lplplotd -lpthread \
	-lfftw3f -llapack -lcblas -latlas \
	-lm"

