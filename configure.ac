#                     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([Camgen], [1.0], [Gijs van den Oord <G.vandenOord@gmail.com>])
AC_CONFIG_SRCDIR([src/unused.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])

# Checks for programs.

OLD_CXXFLAGS=$CXXFLAGS

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL

CXXFLAGS=$OLD_CXXFLAGS

# Automake prequisites

AM_INIT_AUTOMAKE([1.9 -Wall])

# use the C++ compiler for the following checks

AC_LANG([C++])

# Checks for standard header files

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

# Checks for typedefs, structures, and compiler characteristics.

AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for standard library functions.

AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor pow sqrt])

# Debug support (taken from http://www.bioinf.uni-freiburg.de/~mmann/HowTo/automake.html)

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (default=no)])],
    [withdbg="$enableval"],
    [withdbg=no])
AC_MSG_RESULT([$withdbg])

if test "x$withdbg" = "xyes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -Werror -Wno-uninitialized -O0"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O3"
fi

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
if test "x$WITHGNUPLOT" = "xyes"; then
AC_CHECK_PROG(HAVE_GNUPLOT,gnuplot,"yes","no","$GNUPLOT_DIR")
if test "x$HAVE_GNUPLOT" = "xyes"; then
AC_DEFINE_UNQUOTED([GNUPLOTPATH],["$GNUPLOT_DIR/gnuplot"],[The gnuplot binary location])
else
AC_MSG_NOTICE(gnuplot not found in $GNUPLOT_DIR...configuring without.)
fi
fi
])

CAMGEN_CHECK_GNUPLOT

# Check for the LHAPDF header and library

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
if test "x$WITHLHAPDF" = "xyes"; then
AC_CHECK_HEADER([$LHAPDF_DIR/include/LHAPDF/LHAPDF.h],[WITHLHAPDF="yes"],[WITHLHAPDF="no"])
AC_CHECK_LIB([$LHAPDF_DIR/lib/LHAPDF],initPDFSet,[WITHLHAPDF="yes"],[WITHLHAPDF="no"])
AC_CHECK_LIB([$LHAPDF_DIR/lib/LHAPDF],initPDF,[WITHLHAPDF="yes"],[WITHLHAPDF="no"])
AC_CHECK_LIB([$LHAPDF_DIR/lib/LHAPDF],getXmin,[WITHLHAPDF="yes"],[WITHLHAPDF="no"])
AC_CHECK_LIB([$LHAPDF_DIR/lib/LHAPDF],getXmax,[WITHLHAPDF="yes"],[WITHLHAPDF="no"])
AC_CHECK_LIB([$LHAPDF_DIR/lib/LHAPDF],xfx,[WITHLHAPDF="yes"],[WITHLHAPDF="no"])
fi
AC_MSG_RESULT($WITHLHAPDF)
])

CAMGEN_CHECK_LHAPDF

if test "x$WITHLHAPDF"="xyes"; then
    AM_CXXFLAGS="$AM_CXXFLAGS -I${LHAPDF_DIR}/include";
    AM_LDFLAGS="$AM_LDFLAGS -L${LHAPDF_DIR}/lib";
    LIBS = "$LIBS -lLHAPDF";
fi

# Check for the ROOT headers. TODO: check for library functions using custom test progs.

AC_DEFUN([CAMGEN_CHECK_ROOT],[
WITHROOT="yes"
ROOTSYS_DIR="/usr/local"
AC_MSG_CHECKING([for ROOT])
AC_ARG_WITH([ROOT],
	    AC_HELP_STRING([--with-ROOT=DIR],[location of ROOT installation (/usr/local by default).]),
	    [if test "x$with_ROOT" = "xno"; then
	     WITHROOT="no"
	     elif test "x$with_ROOT" != "xyes"; then
	     ROOTSYS_DIR="$with_ROOT";
	     fi],[WITHROOT="no"])
if test "x$WITHROOT" = "xyes"; then
    CPPFLAGS="$CPPFLAGS -I${ROOTSYS_DIR}/include -std=c++11"
    CXXFLAGS="$CXXFLAGS -std=c++11"
    LDFLAGS="$LDFLAGS -L${ROOTSYS_DIR}/lib";
    AC_CHECK_HEADER([TFile.h],[WITHROOT="yes"],[WITHROOT="no"])
    AC_CHECK_HEADER([TTree.h],[WITHROOT="yes"],[WITHROOT="no"])
fi
AC_MSG_RESULT($WITHROOT)
])

CAMGEN_CHECK_ROOT

if test "x$WITHROOT"="xyes"; then
    AM_CXXFLAGS="$AM_CXXFLAGS -I${ROOTSYS_DIR}/include";
    AM_LDFLAGS="$AM_LDFLAGS -L${ROOTSYS_DIR}/lib";
    LIBS="$LIBS -lCore -lTree";
fi

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# Makefiles to configure.

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 lib/Makefile
                 src/Makefile
		 examples/Makefile
                 test/Makefile])
AC_OUTPUT

