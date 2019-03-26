# fscking automake barfs unless a .so is named lib* so this all has to be
# done BY HAND

AM_CPPFLAGS = -I../include

EXTRA_DIST = radiometer.swig

BUILT_SOURCES = radiometer_wrap.c

pkgpython_PYTHON = __init__.py ioutils.py plotutils.py radiometer.py
pkgpyexec_LTLIBRARIES = _radiometer.la

_radiometer_la_SOURCES = radiometer_wrap.c
_radiometer_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS) $(AX_SWIG_PYTHON_CPPFLAGS)
_radiometer_la_CFLAGS = $(AM_CFLAGS) $(PYTHON_CFLAGS) $(AX_SWIG_PYTHON_CFLAGS) $(TCL_CFLAGS)
_radiometer_la_LDFLAGS = -module $(GSL_LIBS) $(FFTW_LIBS) $(TCL_LIBS)
_radiometer_la_LIBADD = ../lib/libradiometer.la

radiometer.py radiometer_wrap.c : radiometer.swig ../include/radiometer/*.h
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I$(SWIG_LIB)/python $(AM_CPPFLAGS) $<
