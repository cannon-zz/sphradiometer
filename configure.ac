AC_INIT([sphradiometer],[0.4.0],[kipp@resceu.s.u-tokyo.ac.jp])
AC_COPYRIGHT([Copyright (C) 2009 Kipp C. Cannon])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AM_INIT_AUTOMAKE([1.9 tar-ustar foreign])
AUTOMAKE_OPTIONS="${AUTOMAKE_OPTIONS} --add-missing --copy"
AC_CONFIG_FILES([ \
	Makefile \
	debian/control \
	debian/Makefile \
	src/Makefile \
	src/include/Makefile \
	src/include/sphradiometer/Makefile \
	src/lib/Makefile \
	src/lib/sphradiometer.pc \
	src/python/Makefile \
	src/bin/Makefile \
	test/Makefile
])

# set pkg-config module path
AC_SUBST([pkgconfigdir],[${libdir}/pkgconfig])

# Checks for programs
AC_PROG_CC_C99
AC_PROG_CXX
AC_OPENMP
AC_PROG_INSTALL
LT_INIT
PKG_PROG_PKG_CONFIG()

# Check for Python
AC_SUBST([MIN_PYTHON_VERSION], [3.0])
AM_PATH_PYTHON([${MIN_PYTHON_VERSION}])
AX_PYTHON_DEVEL()

AC_SUBST([MIN_SWIG_VERSION], [3.0])
AX_PKG_SWIG([${MIN_SWIG_VERSION}],[],[AC_MSG_ERROR("SWIG > ${MIN_SWIG_VERSION} required.")])
AX_SWIG_PYTHON

# Check for math library
AC_CHECK_LIB([m], [main], , [AC_MSG_ERROR([cannot find the math library])])

# Check for FFTW
PKG_CHECK_MODULES([FFTW], [fftw3], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([FFTW_CFLAGS])
AC_SUBST([FFTW_LIBS])

# Check for GSL
AC_SUBST([MIN_GSL_VERSION], [1.9])
PKG_CHECK_MODULES([GSL], [gsl >= ${MIN_GSL_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])

# Check for HDF5
AC_SUBST([MIN_HDF5_VERSION], [1.10])
AX_LIB_HDF5([serial])

# Check for healpix
AC_SUBST([MIN_CHEALPIX_VERSION], [3.30])
PKG_CHECK_MODULES([CHEALPIX], [chealpix], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([CHEALPIX_CFLAGS])
AC_SUBST([CHEALPIX_LIBS])

AC_SUBST([MIN_HEALPIXCXX_VERSION], [3.30])
PKG_CHECK_MODULES([HEALPIXCXX], [healpix_cxx], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([HEALPIXCXX_CFLAGS])
AC_SUBST([HEALPIXCXX_LIBS])

# Check for LAPACK
AC_SUBST([MIN_LAPACK_VERSION], [3.0])
PKG_CHECK_MODULES([LAPACK], [lapack-atlas >= ${MIN_LAPACK_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([LAPACK_CFLAGS])
AC_SUBST([LAPACK_LIBS])
old_CFLAGS="$CFLAGS" CFLAGS="$CFLAGS $LAPACK_CFLAGS"
AC_CHECK_HEADER([lapacke.h], , [AC_MSG_ERROR([cannot find liblapack])])
CFLAGS="$old_CFLAGS"
AC_CHECK_LIB([lapacke], [LAPACKE_zgesdd], , [AC_MSG_ERROR([cannot find liblapack])])

# Check for LAL
AC_SUBST([MIN_LAL_VERSION], [6.19.0])
PKG_CHECK_MODULES([LAL], [lal >= ${MIN_LAL_VERSION} lalsupport >= ${MIN_LAL_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([LAL_CFLAGS])
AC_SUBST([LAL_LIBS])
AC_SUBST([MIN_LALFRAME_VERSION], [1.4.0])
PKG_CHECK_MODULES([LALFRAME], [lalframe >= ${MIN_LALFRAME_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([LALFRAME_CFLAGS])
AC_SUBST([LALFRAME_LIBS])
AC_SUBST([MIN_LALSIMULATION_VERSION], [1.8.0])
PKG_CHECK_MODULES([LALSIMULATION], [lalsimulation >= ${MIN_LALSIMULATION_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([LALSIMULATION_CFLAGS])
AC_SUBST([LALSIMULATION_LIBS])

# Check for S2kit
AC_SUBST([MIN_S2KIT_VERSION], [1.0])
PKG_CHECK_MODULES([S2KIT], [s2kit >= ${MIN_S2KIT_VERSION}], [
	AC_SUBST([HAVE_S2KIT], [1])
], [
	AC_MSG_WARN([Not Found!])
	AC_SUBST([HAVE_S2KIT], [])
])
AM_CONDITIONAL([S2KIT], [test -n $HAVE_S2KIT])
AC_SUBST([S2KIT_CFLAGS])
AC_SUBST([S2KIT_LIBS])

# Check for healpy, numpy, scipy, matplotlib, ohmy
AC_SUBST([MIN_BASEMAP_VERSION], [1.2.0])
AC_SUBST([MIN_HEALPY_VERSION], [1.12.8])
AC_SUBST([MIN_MATPLOTLIB_VERSION], [3.0.0])
AC_SUBST([MIN_NUMPY_VERSION], [1.16.0])
AC_SUBST([MIN_SCIPY_VERSION], [1.1.0])
AX_PYTHON_MODULE([healpy])
AX_PYTHON_MODULE([matplotlib])
AX_PYTHON_MODULE([mpl_toolkits.basemap])
AX_PYTHON_MODULE([numpy], fatal)
AX_PYTHON_MODULE([scipy])
NUMPY_CFLAGS=-I`$PYTHON -c "import numpy;print (numpy.get_include());"`
old_CFLAGS="$CFLAGS" CFLAGS="$CFLAGS ${PYTHON_CPPFLAGS} ${NUMPY_CFLAGS}"
AC_CHECK_HEADER([numpy/arrayobject.h],
	[PYTHON_CPPFLAGS="${PYTHON_CPPFLAGS} ${NUMPY_CFLAGS}"],
	[AC_MSG_ERROR([Numpy extension header not found])],
	[#include "Python.h"])
CFLAGS="$old_CFLAGS"

# CFLAGS, CXXFLAGS
AX_CFLAGS_WARN_ALL([AM_CFLAGS])
AM_CFLAGS="$AM_CFLAGS -Wextra -Wno-missing-field-initializers -Wno-unused-parameter"    # extra gcc-specific stuff
AC_SUBST([AM_CFLAGS])
AX_CXXFLAGS_WARN_ALL([AM_CXXFLAGS])
AM_CXXFLAGS="$AM_CXXFLAGS -Wextra"
AC_SUBST([AM_CXXFLAGS])

# Generate configured output files
AC_PROG_MAKE_SET
AC_OUTPUT
