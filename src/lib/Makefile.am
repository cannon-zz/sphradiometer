AM_CPPFLAGS = -I../include

pkgconfig_DATA = sphradiometer.pc

lib_LTLIBRARIES = libsphradiometer.la libsphinspiralskymap.la

noinst_LTLIBRARIES = libsphradiometerdiag.la

libsphinspiralskymap_la_SOURCES = \
	inspiral_sky_map.c

libsphradiometer_la_SOURCES = \
	correlator.c \
	deconvolution.c \
	inject.c \
	instrument.c \
	projection.c \
	sh_series.c \
	sh_series_deriv.c \
	sh_series_array.c \
	sh_series_healpix.cc \
	sh_series_io.c \
	sh_series_product.c \
	sh_series_rotation.c \
	sh_series_trans.c \
	sky.c

libsphradiometerdiag_la_SOURCES = \
	diagnostics.c

libsphinspiralskymap_la_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS) $(LAPACK_CFLAGS)
libsphinspiralskymap_la_CXXFLAGS = $(AM_CXXFLAGS)
libsphinspiralskymap_la_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS) $(LAPACK_LIBS)
libsphinspiralskymap_la_LIBADD = libsphradiometer.la

libsphradiometer_la_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS) $(LAPACK_CFLAGS)
libsphradiometer_la_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CFLAGS) $(HEALPIXCXX_CFLAGS)
libsphradiometer_la_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS) $(LAPACK_LIBS) $(HEALPIXCXX_LIBS)

libsphradiometerdiag_la_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) $(LALSIMULATION_CFLAGS) $(LAPACK_CFLAGS)
libsphradiometerdiag_la_CXXFLAGS = $(AM_CXXFLAGS)
libsphradiometerdiag_la_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAL_LIBS) $(LALSIMULATION_LIBS) $(LAPACK_LIBS) $(HEALPIXCXX_LIBS)
libsphradiometerdiag_la_LIBADD = libsphradiometer.la
