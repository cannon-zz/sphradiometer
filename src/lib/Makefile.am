INCLUDES = -I../include

pkgconfig_DATA = radiometer.pc

lib_LTLIBRARIES = libradiometer.la

libradiometer_la_SOURCES = sky.c inject.c sh_series.c sh_series_array.c sh_series_product.c sh_series_rotation.c instrument.c projection.c correlator.c deconvolution.c
libradiometer_la_LDFLAGS = $(FFTW_LIBS) $(GSL_LIBS)
