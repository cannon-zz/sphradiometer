INCLUDES = -I../include -I.

bin_PROGRAMS = radiometer fastX
noinst_PROGRAMS = tests

radiometer_SOURCES = main.c instruments.c diagnostics.c diagnostics.h instruments.h
radiometer_LDFLAGS = $(GSL_LIBS) $(LAL_LIBS)
radiometer_LDADD = ../lib/libradiometer.la

fastX_SOURCES = fastX.c instruments.c instruments.h
fastX_LDFLAGS = $(GSL_LIBS) $(LAL_LIBS) $(LALFRAME_LIBS)
fastX_LDADD = ../lib/libradiometer.la

tests_SOURCES = tests.c instruments.c diagnostics.c diagnostics.h instruments.h
tests_LDFLAGS = $(GSL_LIBS) $(LAL_LIBS)
tests_LDADD = ../lib/libradiometer.la
