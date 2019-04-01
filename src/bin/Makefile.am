AM_CPPFLAGS = -I../include -I.

bin_PROGRAMS = inspiral_sky_map radiometer fastX
noinst_PROGRAMS = tests

radiometer_SOURCES = main.c instruments.c instruments.h
radiometer_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
radiometer_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
radiometer_LDADD = ../lib/libsphradiometer.la

fastX_SOURCES = fastX.c instruments.c instruments.h
fastX_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALFRAME_CFLAGS) $(LALSIMULATION_CFLAGS)
fastX_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALFRAME_LIBS) $(LALSIMULATION_LIBS)
fastX_LDADD = ../lib/libsphradiometer.la

inspiral_sky_map_SOURCES = inspiral_sky_map.c instruments.c instruments.h
inspiral_sky_map_CFLAGS = $(AM_CFLAGS) $(CHEALPIX_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS) $(LAL_CFLAGS) $(LALFRAME_CFLAGS) $(LALSIMULATION_CFLAGS)
inspiral_sky_map_LDFLAGS = $(AM_LDFLAGS) $(CHEALPIX_LIBS) $(FFTW_LIBS) $(GSL_LIBS) $(HDF5_LIBS) $(LAL_LIBS) $(LALFRAME_LIBS) $(LALSIMULATION_LIBS)
inspiral_sky_map_LDADD = ../lib/libsphradiometer.la

tests_SOURCES = tests.c
tests_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
tests_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
tests_LDADD = ../lib/libsphradiometer.la
