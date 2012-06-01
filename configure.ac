AC_INIT([rings], [1.0], [w-farr@northwestern.edu])
AM_INIT_AUTOMAKE

AC_PROG_CC
LT_INIT()

AX_PATH_GSL(1.12)

AC_CHECK_FUNC([getopt_long], [], [AC_MSG_ERROR([could not call getopt_long])])

AC_CHECK_HEADER([getopt.h], [], [AC_MSG_ERROR([could not find getopt.h])])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([
        Makefile
        src/Makefile
        test/Makefile
	doc/Makefile
])

AC_OUTPUT