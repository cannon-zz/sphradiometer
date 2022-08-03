AM_CPPFLAGS = -I../include -I.

bin_PROGRAMS = inspiral_sky_map fastX

noinst_PROGRAMS = radiometer

radiometer_SOURCES = main.c
radiometer_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS)
radiometer_LDFLAGS = $(AM_LDFLAGS) $(LAL_LIBS) $(FFTW_LIBS) $(GSL_LIBS)
radiometer_LDADD = ../lib/libsphradiometer.la ../lib/libsphradiometerdiag.la

fastX_SOURCES = fastX.c
fastX_CFLAGS = $(AM_CFLAGS) $(LALSIMULATION_CFLAGS) $(LALFRAME_CFLAGS) $(LAL_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS)
fastX_LDFLAGS = $(AM_LDFLAGS) $(LALSIMULATION_LIBS) $(LALFRAME_LIBS) $(LAL_LIBS) $(FFTW_LIBS) $(GSL_LIBS)
fastX_LDADD = ../lib/libsphradiometer.la

inspiral_sky_map_SOURCES = inspiral_sky_map.c
inspiral_sky_map_CFLAGS = $(AM_CFLAGS) $(LALSIMULATION_CFLAGS) $(LALFRAME_CFLAGS) $(LAL_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS)
inspiral_sky_map_LDFLAGS = $(AM_LDFLAGS) $(LALSIMULATION_LIBS) $(LALFRAME_LIBS) $(LAL_LIBS) $(FFTW_LIBS) $(GSL_LIBS) $(HDF5_LDFLAGS) $(HDF5_LIBS)
inspiral_sky_map_LDADD = ../lib/libsphinspiralskymap.la ../lib/libsphradiometer.la
