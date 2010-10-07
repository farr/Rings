AC_INIT([rings], [1.0], [w-farr@northwestern.edu])
AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_LIBTOOL

AC_ARG_WITH([gsl], 
            AS_HELP_STRING([--with-gsl=DIR], [directory of the GNU Scientific Library]),
            [CFLAGS="-I$withval/include $CFLAGS"; CPPFLAGS="-I$withval/include $CPPFLAGS"; LDFLAGS="-L$withval/lib $LDFLAGS"])

AC_CHECK_FUNC([getopt_long], [], [AC_MSG_ERROR([could not call getopt_long])])

AC_CHECK_LIB([gslcblas], [cblas_ddot], 
             [AC_DEFINE([HAVE_LIBGSLCBLAS], [], [GSL CBLAS]) LIBS="-lgslcblas $LIBS"],
             [AC_MSG_ERROR([could not find GSL CBLAS library])])
AC_CHECK_LIB([gsl], [gsl_sf_ellint_Ecomp],
             [AC_DEFINE([HAVE_LIBGSL], [], [GSL]) LIBS="-lgsl $LIBS"],
             [AC_MSG_ERROR([could not find GSL library])])

AC_CHECK_HEADER([gsl/gsl_sf.h], [], [AC_MSG_ERROR([could not find GSL headers; consider using --with-gsl argument])])
AC_CHECK_HEADER([getopt.h], [], [AC_MSG_ERROR([could not find getopt.h])])

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
        Makefile
        src/Makefile
        test/Makefile
])

AC_OUTPUT