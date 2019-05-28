AM_CPPFLAGS = -I../include

pkgconfig_DATA = sphradiometer.pc

lib_LTLIBRARIES = libsphradiometer.la

libsphradiometer_la_SOURCES = \
	correlator.c \
	deconvolution.c \
	diagnostics.c \
	inject.c \
	instrument.c \
	projection.c \
	sh_series.c \
	sh_series_array.c \
	sh_series_healpix.cc \
	sh_series_io.c \
	sh_series_product.c \
	sh_series_rotation.c \
	sh_series_trans.c \
	sky.c

libsphradiometer_la_CFLAGS = $(AM_CFLAGS) $(CHEALPIX_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS) $(LAPACK_CFLAGS)
libsphradiometer_la_CXXFLAGS = $(AM_CXXFLAGS) $(HEALPIXCXX_CFLAGS)
libsphradiometer_la_LDFLAGS = $(AM_LDFLAGS) $(CHEALPIX_LIBS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS) $(LAPACK_LIBS) $(HEALPIXCXX_LIBS)
