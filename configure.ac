# configure.in for PSFEx
# (C) E.Bertin 2002-2007
# Process this file with autoconf to produce a configure script.
# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT(psfex, 2.0.9, [bertin@iap.fr])
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_AUX_DIR(autoconf)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "Emmanuel Bertin")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_atlas.m4)
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_urbi_resolve_dir.m4)

# Display pakage and version number
AC_MSG_RESULT([*********** Configuring:  $PACKAGE $VERSION ($date) **********])

# Initialize the list of compilers to consider
cclist="cc gcc"

# Backup and reset the input CFLAGS and LDFLAGS
mycflags="$CFLAGS"
CFLAGS=""
myldflags="$LDFLAGS"
LDFLAGS=""

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([for Linux Intel C compiler mode])
AC_ARG_ENABLE(icc,
	[AC_HELP_STRING([--enable-icc],
	[Enable special mode for compilation with the Intel compiler \
(off by default)])],
        use_icc="yes"
        cclist="icc $cclist"
	AC_MSG_RESULT([yes]),
        use_icc="no"
	AC_MSG_RESULT([no]))

# Checks for programs.
# GCC is chosen last because it is likely to yield less optimized code
AC_LANG(C)
AC_PROG_CC([$cclist])
# C Compiler: Check that it is ANSI C
AM_PROG_CC_STDC
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

# C Compiler: Check that it is POSIX-compliant
AC_ISC_POSIX
ACX_PROG_CC_OPTIM
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stdlib.h string.h sys/mman.h sys/types.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv memcpy memmove memset mkdir munmap strstr \
sysconf])

# Check support for large files
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Set the data directory to a true absolute path
datadir2=$(URBI_RESOLVE_DIR([$datadir]))

# Provide a special option for ATLAS
AC_ARG_WITH(atlas_libpath,
	[AC_HELP_STRING([--with-atlas_libpath=<ATLAS library path>],
	[Provide an alternative path to the ATLAS library])],
    atlas_libdir=$withval,
    atlas_libdir=""
    )

# Provide a special option for the default XSLT URL
AC_ARG_WITH(xsl_url,
        [AC_HELP_STRING([--with-xsl_url=<default URL for XSLT filter>],
        [Provide an alternative default URL of the XSLT filter])],
    [xsl_url=$withval],
    [xsl_url="file://$datadir2/$PACKAGE_NAME/$PACKAGE_NAME.xsl"]
    )

AC_DEFINE_UNQUOTED([XSL_URL], "$xsl_url",[Default URL of the XSLT filter])

# Provide special option for gprof profiling
AC_MSG_CHECKING([for gprof profiler mode])
AC_ARG_ENABLE(gprof,
	[AC_HELP_STRING([--enable-gprof],
	[Enable special mode for compilation with the gprof profiler \
(off by default)])],
        use_gprof="yes"
	AC_MSG_RESULT([yes]),
        use_gprof="no"
	AC_MSG_RESULT([no]))

# Enable static linking
AC_MSG_CHECKING([for static linking])
AC_ARG_ENABLE(static,
	[AC_HELP_STRING([--enable-static],
	[Enable static linking (off by default)])],
        use_static="$enableval"
	AC_MSG_RESULT([yes]),
        use_static="no"
	AC_MSG_RESULT([no]))

ACX_ATLAS($atlas_libdir,[use_atlas=yes],[use_atlas=no])

################## handle the ATLAS library(linear algebra) ##################
if test "$use_atlas" = "yes"; then
  if test "$use_pthreads" = "yes"; then
    LIBS="-llapack -lptcblas -lcblas -latlas $LIBS"
  else
    LIBS="-llapack -lcblas -latlas $LIBS"
  fi
else
  AC_MSG_ERROR([CBLAS/LAPack library not found! Exiting.])
fi

# Special optimization options for the INTEL C compiler
if test $use_icc = "yes" -a $CC = "icc"; then
  CFLAGS="-O -axWP -ip -unroll"
fi

if test "$use_gprof" = "yes"; then
  if test "$use_icc" = "yes"; then
    CFLAGS="$CFLAGS -pq"
  else
    CFLAGS="$CFLAGS -pg"
  fi
  use_static="no"
fi

if test "$use_static" = "yes"; then
  LDFLAGS="$LDFLAGS -all-static"
fi

# Override automatic CFLAGS and LDFLAGS with those of user
if test -n "$mycflags"; then
CFLAGS="$mycflags"
fi
if test -n "$myldflags"; then
LDFLAGS="$myldflags"
fi

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
man/Makefile
src/Makefile
src/fits/Makefile
src/levmar/Makefile
psfex.spec
man/psfex.1])
AC_OUTPUT
