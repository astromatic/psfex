dnl @synopsis ACX_PLPLOT([PLPLOT_DIR,[ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the PlPlot library and header files
dnl are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$PLPLOT_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $PLPLOT_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if PlPlot
dnl is found (HAVE_PLPLOT is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.
dnl
dnl @version $Id: acx_plplot.m4,v 1.0 2004/05/30 21:30:17 bertin Exp $
dnl @author Emmanuel Bertin <bertin@iap.fr>

AC_DEFUN([ACX_PLPLOT], [
AC_REQUIRE([AC_CANONICAL_HOST])

AC_CHECK_PROG(acx_plplot_ok, plplot-config, [yes], [no])
plpath=`plplot-config --prefix`
if test x$acx_plplot_ok = xyes; then
AC_CHECK_HEADER([${plpath}/include/plplot/plplot.h],[acx_plplot_ok=yes])
fi
if test x$acx_plplot_ok = xyes; then
    [PLPLOT_CFLAGS=`plplot-config --cflags`]
    [PLPLOT_DIR="${plpath}"]
    [PLPLOT_LIBPATH="-L${plpath}/lib"]
    [PLPLOT_LIBS="-lplplotd"]
fi

AC_SUBST(PLPLOT_LIBS)
AC_SUBST(PLPLOT_CFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_plplot_ok" = xyes; then
        AC_DEFINE(HAVE_PLPLOT,1,
        [Define if you have the PLPlot libraries and header files.])
        $2
else
        $3
fi

])dnl ACX_PLPLOT
