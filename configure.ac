AC_INIT([radiometer],[0.3.0],[kipp@gravity.phys.uwm.edu])
AC_COPYRIGHT([Copyright (C) 2009 Kipp C. Cannon])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_FILES([ \
	Makefile \
	src/Makefile \
	src/include/Makefile \
	src/include/radiometer/Makefile \
	src/lib/Makefile \
	src/lib/radiometer.pc \
	src/python/Makefile \
	src/bin/Makefile \
	src/test/Makefile \
	src/matlab/Makefile \
])

# Checks for programs
AC_PROG_CC_C99
AC_PROG_INSTALL
LT_INIT
PKG_PROG_PKG_CONFIG()

# Check for Python
AM_PATH_PYTHON([2.3])
AX_PKG_SWIG([1.3.24])
SWIG_PYTHON

# Check for math library
AC_CHECK_LIB([m], [main], , [AC_MSG_ERROR([cannot find the math library])])

# Check for FFTW
PKG_CHECK_MODULES([FFTW], [fftw3], , [echo "Not Found!" ; exit 1])
AC_SUBST([FFTW_CFLAGS])
AC_SUBST([FFTW_LIBS])

# Check for GSL
PKG_CHECK_MODULES([GSL], [gsl], , [echo "Not Found!" ; exit 1])
AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])

# Check for LAL
PKG_CHECK_MODULES([LAL], [lal], , [echo "Not Found!" ; exit 1])
AC_SUBST([LAL_CFLAGS])
AC_SUBST([LAL_LIBS])
PKG_CHECK_MODULES([LALFRAME], [lalframe], , [echo "Not Found!" ; exit 1])
AC_SUBST([LALFRAME_CFLAGS])
AC_SUBST([LALFRAME_LIBS])
PKG_CHECK_MODULES([LALSIMULATION], [lalsimulation], , [echo "Not Found!" ; exit 1])
AC_SUBST([LALSIMULATION_CFLAGS])
AC_SUBST([LALSIMULATION_LIBS])

# pkg-config module path
AC_SUBST([pkgconfigdir],[${libdir}/pkgconfig])

# Set package-wide CFLAGS
AC_SUBST([AM_CFLAGS], ["-std=c99 -Wall -Wno-unused-function \$(FFTW_CFLAGS) \$(GSL_CFLAGS) \$(LAL_CFLAGS)"])

# Generate configured output files
AC_PROG_MAKE_SET
AC_OUTPUT
