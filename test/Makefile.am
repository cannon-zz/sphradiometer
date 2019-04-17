AM_CPPFLAGS = -I../src/include

TESTS = \
	correlator_validate \
	healpy_validate.py \
	s2kit_validate \
	sh_series_validate

check_PROGRAMS = \
	correlator_validate \
	s2kit_validate \
	sh_series_validate

correlator_validate_SOURCES = correlator_validate.c
correlator_validate_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
correlator_validate_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
correlator_validate_LDADD = ../src/lib/libsphradiometer.la

s2kit_validate_SOURCES = s2kit_validate.c
s2kit_validate_CFLAGS = $(AM_CFLAGS) $(S2KIT_CFLAGS)
s2kit_validate_LDFLAGS = $(AM_LDLAGS) $(S2KIT_LIBS)
s2kit_validate_LDADD = ../src/lib/libsphradiometer.la

sh_series_validate_SOURCES = sh_series_validate.c
sh_series_validate_CFLAGS = $(AM_CFLAGS)
sh_series_validate_LDFLAGS = $(AM_LDLAGS)
sh_series_validate_LDADD = ../src/lib/libsphradiometer.la

clean-local :
	-rm -f tests_{exact,fdaverage,tdaverage}_*.dat
