AC_INIT([sphradiometer],[0.5.0],[kipp@resceu.s.u-tokyo.ac.jp])
AC_COPYRIGHT([Copyright (C) 2009 Kipp C. Cannon])
AC_CONFIG_HEADERS([src/include/sphradiometer/sphradiometer_config.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AM_INIT_AUTOMAKE([1.9 tar-ustar foreign])
AUTOMAKE_OPTIONS="${AUTOMAKE_OPTIONS} --add-missing --copy"
AC_CONFIG_FILES([ \
	Makefile \
	sphradiometer.spec \
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
# FIXME:  settle on AC_PROG_CC when LDG cluters have autoconf >= 2.70
m4_version_prereq([2.70], [AC_PROG_CC()], [AC_PROG_CC_C99])
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
# minimum required version is 2.0 (the version that added gsl_spline2d.h
# and switched to native C complex types)
AC_SUBST([MIN_GSL_VERSION], [2.0])
PKG_CHECK_MODULES([GSL], [gsl >= ${MIN_GSL_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])

# Check for HDF5
AC_SUBST([MIN_HDF5_VERSION], [1.8])
AX_LIB_HDF5([serial])

# Check for healpix
AC_SUBST([MIN_HEALPIXCXX_VERSION], [3.30])
PKG_CHECK_MODULES([HEALPIXCXX], [healpix_cxx], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([HEALPIXCXX_CFLAGS])
AC_SUBST([HEALPIXCXX_LIBS])

# Check for LAPACK.  the following has been tested to work on Debian. (and
# also on ldas-pcdev1.ligo.caltech.edu as of 2022-08-02).  YMMV.
AS_IF([test -f /etc/redhat-release], [
	# redhat packages frequently fail to provide pkg-config metadata.
	# I don't have the time needed to fix redhat's problems, so I'm
	# just going to hard code this.  if this doesn't work for you, feel
	# free to fight with this yourself.
	LAPACK_CFLAGS=-"I/usr/include/openblas"
	LAPACK_LIBS="-lopenblaso"
], [
	# on Debian, lapacke is its own package and acts as a front-end for
	# whichever of the BLAS implementations you choose to install.
	# testing seems to show that if lapacke-dev is installed, and we
	# ask pkg-config for "lapack" we get CFLAGS and LIBS that pick up
	# the lapacke.h from lapacke-dev and gives you whichever BLAS
	# library you have installed (I have one of the openblas variants).
	PKG_CHECK_MODULES([LAPACK], [lapack], , [AC_MSG_ERROR([Not Found!])])
])
AC_SUBST([LAPACK_CFLAGS])
AC_SUBST([LAPACK_LIBS])
ac_save_CFLAGS=$CFLAGS
CFLAGS=$LAPACK_CFLAGS
AC_CHECK_HEADER([lapacke.h], , [AC_MSG_ERROR([cannot find lapacke.h])])
CFLAGS=$ac_save_CFLAGS

# Check for LAL
AC_SUBST([MIN_LAL_VERSION], [7.0.0])
PKG_CHECK_MODULES([LAL], [lal >= ${MIN_LAL_VERSION} lalsupport >= ${MIN_LAL_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([LAL_CFLAGS])
AC_SUBST([LAL_LIBS])
AC_SUBST([MIN_LALFRAME_VERSION], [2.0.0])
PKG_CHECK_MODULES([LALFRAME], [lalframe >= ${MIN_LALFRAME_VERSION}], , [AC_MSG_ERROR([Not Found!])])
AC_SUBST([LALFRAME_CFLAGS])
AC_SUBST([LALFRAME_LIBS])
AC_SUBST([MIN_LALSIMULATION_VERSION], [3.0.0])
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
