#                     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([Camgen], [1.0], [Gijs van den Oord <G.vandenOord@gmail.com>])
AC_CONFIG_SRCDIR([src/unused.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])

# Checks for programs.

AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Automake prequisites

AM_INIT_AUTOMAKE([1.9])

# use the C++ compiler for the following checks

AC_LANG([C++])

# Checks for header files

AC_HEADER_STDC
AC_CHECK_HEADER([string])
AC_CHECK_HEADER([iostream])
AC_CHECK_HEADER([fstream])
AC_CHECK_HEADER([sstream])
AC_CHECK_HEADER([bitset])
AC_CHECK_HEADER([list])
AC_CHECK_HEADER([map])
AC_CHECK_HEADER([set])
AC_CHECK_HEADER([vector])
AC_CHECK_HEADER([algorithm])

# Check for the LHAPDF header.

AC_DEFUN([CAMGEN_CHECK_LHAPDF],[
WITHLHAPDF="yes"
LHAPDF_DIR="/usr/local"
AC_MSG_CHECKING([for LHAPDF])
AC_ARG_WITH([LHAPDF],
	    AC_HELP_STRING([--with-LHAPDF=DIR],[location of LHAPDF installation (/usr/local by default).]),
	    [if test "x$with_LHAPDF" = "xno"; then
	     WITHLHAPDF="no"
	     elif test "x$with_LHAPDF" != "xyes"; then
	     LHAPDF_DIR="$with_LHAPDF";
	     fi],[WITHLHAPDF="no"])
AC_MSG_RESULT($WITHLHAPDF)
if test "$WITHLHAPDF" = "yes"; then
LDFLAGS="-L${LHAPDF_DIR}/lib $LDFLAGS";
CPPFLAGS="-I${LHAPDF_DIR}/include $CPPFLAGS";
AC_CHECK_HEADER([$LHAPDF_DIR/include/LHAPDF/LHAPDF.h],
		 AC_DEFINE([HAVE_LHAPDF_H_],[1],[Define 1 if you have the LHAPDF header installed.]),
		 AC_MSG_NOTICE(LHAPDF.h not found in $LHAPDF_DIR/include/LHAPDF/ 
			       or not compilable with flags $CPPFLAGS...configuring without.),[])
fi
])
CAMGEN_CHECK_LHAPDF
AM_CONDITIONAL([USE_LHAPDF],[test "$WITHLHAPDF" = "yes"])
AC_SUBST(LHAPDF_DIR)

# Checks for typedefs, structures, and compiler characteristics.

AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.

AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor pow sqrt])

# Check for gnuplot

AC_DEFUN([CAMGEN_CHECK_GNUPLOT],[
WITHGNUPLOT="yes"
GNUPLOT_DIR="/usr/local/bin"
AC_ARG_WITH([gnuplot],
	    AC_HELP_STRING([--with-gnuplot=DIR],[location of gnuplot
			    installation (/usr/local/bin by default)]),
	    [if test "x$with_gnuplot" = "xno"; then
	     WITHGNUPLOT="no"
	     elif test "x$with_gnuplot" != "xyes"; then
	     GNUPLOT_DIR="$with_gnuplot"
	     fi],[WITHGNUPLOT="no"])
AC_MSG_CHECKING([for gnuplot location])
AC_MSG_RESULT($WITHGNUPLOT)
if test "$WITHGNUPLOT" = "yes"; then
AC_CHECK_PROG(HAVE_GNUPLOT,gnuplot,"yes","no","$GNUPLOT_DIR")
if test "$HAVE_GNUPLOT" = "yes"; then
AC_DEFINE_UNQUOTED([GNUPLOTPATH],["$GNUPLOT_DIR/gnuplot"],[The gnuplot binary location])
else
AC_MSG_NOTICE(gnuplot not found in $GNUPLOT_DIR...configuring without.)
fi
fi
])
CAMGEN_CHECK_GNUPLOT

# Makefiles to configure.

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 lib/Makefile
                 src/Makefile
		 examples/Makefile
                 test/Makefile])
AC_OUTPUT

