AM_CPPFLAGS = -I../src/include

TESTS = \
	correlator_validate \
	healpy_validate.py \
	rotation_validate \
	s2kit_validate \
	sh_series_validate

check_PROGRAMS = \
	correlator_validate \
	rotation_validate \
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

rotation_validate_SOURCES = rotation_validate.c
rotation_validate_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
#rotation_validate_LDFLAGS = $(AM_LDFLAGS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
rotation_validate_LDADD = ../src/lib/libsphradiometer.la

sh_series_validate_SOURCES = sh_series_validate.c
sh_series_validate_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
#sh_series_validate_LDFLAGS = $(AM_LDFLAGS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
sh_series_validate_LDADD = ../src/lib/libsphradiometer.la

plots :
	-rm -vf tests_*.png
	./plottests tests_*.fits

clean-local :
	-rm -f tests_{exact,fdaverage,tdaverage}_*.fits
	-rm -f tests_*.png
