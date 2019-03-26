AM_CPPFLAGS = -I../include -I.

bin_PROGRAMS = radiometer fastX
noinst_PROGRAMS = tests inspiral_sky_map

radiometer_SOURCES = main.c instruments.c diagnostics.c diagnostics.h instruments.h
radiometer_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
radiometer_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
radiometer_LDADD = ../lib/libsphradiometer.la

fastX_SOURCES = fastX.c instruments.c instruments.h
fastX_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALFRAME_CFLAGS) $(LALSIMULATION_CFLAGS)
fastX_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS) $(LAL_LIBS) $(LALFRAME_LIBS) $(LALSIMULATION_LIBS)
fastX_LDADD = ../lib/libsphradiometer.la

inspiral_sky_map_SOURCES = inspiral_sky_map.c instruments.c instruments.h diagnostics.c diagnostics.h
inspiral_sky_map_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALFRAME_CFLAGS) $(LALSIMULATION_CFLAGS)
inspiral_sky_map_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS) $(LAL_LIBS) $(LALFRAME_LIBS) $(LALSIMULATION_LIBS)
inspiral_sky_map_LDADD = ../lib/libsphradiometer.la

tests_SOURCES = tests.c instruments.c diagnostics.c diagnostics.h instruments.h
tests_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
tests_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
tests_LDADD = ../lib/libsphradiometer.la
