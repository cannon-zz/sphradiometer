AM_CPPFLAGS = -I../include

# This is a trick taken from the gst-python automake setup.
# All of the Python scripts will be installed under the exec dir,
# which prevents the module from getting spread across lib and lib64
# on e.g. CentOS.
pkgpythondir = $(pkgpyexecdir)

EXTRA_DIST = sphradiometer.swig

BUILT_SOURCES = sphradiometer.py sphradiometer_wrap.c

pkgpython_PYTHON = __init__.py ioutils.py plotutils.py sphradiometer.py
pkgpyexec_LTLIBRARIES = _sphradiometer.la

_sphradiometer_la_SOURCES = sphradiometer_wrap.c
_sphradiometer_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS) $(AX_SWIG_PYTHON_CPPFLAGS)
_sphradiometer_la_CFLAGS = $(AM_CFLAGS) $(PYTHON_CFLAGS) $(AX_SWIG_PYTHON_CFLAGS) $(TCL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
_sphradiometer_la_LDFLAGS = -module $(GSL_LIBS) $(FFTW_LIBS) $(TCL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
_sphradiometer_la_LIBADD = ../lib/libsphradiometer.la

sphradiometer.py sphradiometer_wrap.c : sphradiometer.swig ../include/sphradiometer/*.h
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I$(SWIG_LIB)/python $(AM_CPPFLAGS) $<
