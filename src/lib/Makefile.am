AM_CPPFLAGS = -I../include

pkgconfig_DATA = radiometer.pc

lib_LTLIBRARIES = libradiometer.la

libradiometer_la_SOURCES = sky.c inject.c sh_series.c sh_series_array.c sh_series_product.c sh_series_rotation.c instrument.c projection.c correlator.c deconvolution.c
libradiometer_la_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS) $(GSL_CFLAGS) $(LAPACK_CFLAGS)
libradiometer_la_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS) $(GSL_LIBS) $(LAPACK_LIBS)
