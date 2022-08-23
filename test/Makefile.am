AM_CPPFLAGS = -I${srcdir}/../src/include

EXTRA_DIST = \
	healpy_validate.py

TESTS = \
	correlator_validate \
	deriv_validate \
	healpy_validate.py \
	product_validate \
	rotation_validate \
	sh_series_validate \
	sh_series_io_validate

check_PROGRAMS = \
	correlator_validate \
	deriv_validate \
	product_validate \
	rotation_validate \
	sh_series_validate \
	sh_series_io_validate

if HAVE_GSL_2_0
TESTS += \
	gsl_validate
check_PROGRAMS += \
	gsl_validate
endif

if S2KIT
TESTS += \
	s2kit_validate
check_PROGRAMS += \
	s2kit_validate
endif

correlator_validate_SOURCES = correlator_validate.c
correlator_validate_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
correlator_validate_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
correlator_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

deriv_validate_SOURCES = deriv_validate.c
deriv_validate_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
deriv_validate_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
deriv_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

gsl_validate_SOURCES = gsl_validate.c
gsl_validate_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
gsl_validate_LDFLAGS = $(AM_LDLAGS) $(GSL_LIBS)
gsl_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

product_validate_SOURCES = product_validate.c
product_validate_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
#product_validate_LDFLAGS = $(AM_LDFLAGS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
product_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

rotation_validate_SOURCES = rotation_validate.c
rotation_validate_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(GSL_CFLAGS) $(LALSIMULATION_CFLAGS)
rotation_validate_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS)
rotation_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

s2kit_validate_SOURCES = s2kit_validate.c
s2kit_validate_CFLAGS = $(AM_CFLAGS) $(S2KIT_CFLAGS)
s2kit_validate_LDFLAGS = $(AM_LDLAGS) $(S2KIT_LIBS)
s2kit_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

sh_series_validate_SOURCES = sh_series_validate.c
sh_series_validate_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
#sh_series_validate_LDFLAGS = $(AM_LDFLAGS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
sh_series_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

sh_series_io_validate_SOURCES = sh_series_io_validate.c
sh_series_io_validate_CFLAGS = $(AM_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS)
#sh_series_io_validate_LDFLAGS = $(AM_LDFLAGS) $(LAL_LIBS) $(LALSIMULATION_LIBS)
sh_series_io_validate_LDADD = ../src/lib/libsphradiometer.la ../src/lib/libsphradiometerdiag.la

plots :
	-rm -vf tests_*.png
	./plottests tests_*.fits

clean-local :
	-rm -f test.fits
	-rm -f tests_{exact,fdaverage,tdaverage}_*.fits
	-rm -f tests_*.png
