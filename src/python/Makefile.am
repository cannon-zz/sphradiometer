# fscking automake barfs unless a .so is named lib* so this all has to be
# done BY HAND

INCLUDES = -I../include

EXTRA_DIST = radiometer.swig

BUILT_SOURCES = radiometer_wrap.c

pkgpython_PYTHON = __init__.py ioutils.py plotutils.py radiometer.py
pkgpyexec_LTLIBRARIES = _radiometer.la

_radiometer_la_SOURCES = radiometer_wrap.c
_radiometer_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS)
_radiometer_la_LDFLAGS = -module $(GSL_LIBS) $(FFTW_LIBS)
_radiometer_la_LIBADD = ../lib/libradiometer.la

radiometer.py radiometer_wrap.c : radiometer.swig ../include/radiometer/*.h
	$(SWIG) $(SWIG_PYTHON_OPT) $(INCLUDES) $<
