AM_CPPFLAGS = -I../src/include

TESTS = \
	healpy_validate \
	s2kit_validate \
	tests

check_PROGRAMS = \
	s2kit_validate \
	tests

s2kit_validate_SOURCES = s2kit_validate.c
s2kit_validate_CFLAGS = $(AM_CFLAGS) $(S2KIT_CFLAGS)
s2kit_validate_LDFLAGS = $(AM_LDLAGS) $(S2KIT_LIBS)
s2kit_validate_LDADD = ../src/lib/libsphradiometer.la

tests_SOURCES = tests.c
tests_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
tests_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
tests_LDADD = ../src/lib/libsphradiometer.la

clean-local :
	-rm -f tests_{exact,fdaverage,tdaverage}_*.dat
