# Copyright (C) 2019  Kipp C. Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


import healpy
import math
import numpy


from . import sphradiometer


#
# =============================================================================
#
#                                  Conversion
#
# =============================================================================
#


def healpy_alm_l_max(alm, m_max):
	"""
	Compute l_max from the length of the healpy alm array, given the
	value of m_max.  Raises ValueError if the (l_max, m_max)
	combination is not compatible with the requirements of the
	sphradiometer library (even if healpix allows it).

	NOTE:  healpix does not support complex-valued functions on the
	sphere, so the negative m-coefficients are not present in the alm
	array.  The relationship between the length of the array and the
	maximum spherical harmonic order it contains is not, in general,
	the same as the relationship between these things for the
	sphradiometer library.
	"""
	n, = alm.shape	# confirm 1-D
	l_max = healpy.Alm.getlmax(n, m_max)
	if m_max not in (0, l_max):
		raise ValueError("unsupported (l_max, m_max) = (%d, %d)" % (l_max, m_max))
	return l_max


def healpy_alm_length(l_max, m_max):
	"""
	Compute the length of the healpy alm array for the given l_max and
	m_max parameters.  Raises ValueError if the (l_max, m_max)
	combination is not compatible with the requirements of the
	sphradiometer library (even if healpix allows it).

	NOTE:  healpix does not support complex-valued functions on the
	sphere, so the negative m-coefficients are not present in the alm
	array.  The relationship between the length of the array and the
	maximum spherical harmonic order it contains is not, in general,
	the same as the relationship between these things for the
	sphradiometer library.
	"""
	if m_max not in (0, l_max):
		raise ValueError("unsupported (l_max, m_max) = (%d, %d)" % (l_max, m_max))
	return healpy.Alm.getsize(l_max, m_max)


def sh_series_to_healpy_alm(series):
	"""
	Convert an sh_series object to a healpy-compatible alm array.  The
	return value is a 3-element tuple (alm, l_max, m_max).

	NOTE:  healpix does not support complex-valued functions on the
	sphere, so the negative-m coefficients are not present in the alm
	array.  There is no check performed to confirm that series
	describes a real-valued function.  It is the obligation of the
	calling code to ensure this.  The negative-m coefficients will be
	silently discarded by this function.
	"""
	l_max = series.l_max
	m_max = 0 if series.polar else l_max
	alms = numpy.ndarray(shape = (healpy_alm_length(l_max, m_max),), dtype = "complex128")
	# the coefficient order for healpix and sphradiometer is identical.
	# this can be confirmed by running the library's test suite, which
	# contains a test to confirm this, index-by-index, for each
	# coefficient.  healpix stores only the non-negative m
	# coefficientes, which are at the start of the array, so after
	# determining the length of the healpix array, we copy one-to-one
	# the first that-many coefficients from the sh_series object.
	for i in range(len(alms)):
		alms[i] = sphradiometer.double_complex_array_getitem(series.coeff, i)
	return alms, l_max, m_max


def healpy_alm_to_sh_series(alm, l_max, m_max):
	"""
	Convert a healpy-compatible alm array to an sh_series object.

	NOTE:  healpix does not support complex-valued functions on the
	sphere, so the negative-m coefficients are not present in the alm
	array.
	"""
	# we could simply compute l_max from m_max and the length of alm,
	# and since we want to safety check the input we have to recompute
	# l_max anyway, so it might make more sense to not require it be
	# part of the input.  I think it provides better API consistency if
	# we always deal with the (alm, l_max, m_max) triple.
	assert l_max == healpy_alm_l_max(alm, m_max), "invalid l_max=%d" % l_max

	series = sh_series_new(l_max, m_max == 0)

	for m in range(0, m_max + 1):
		for l in range(m, l_max + 1):
			x = alm[healpy.Alm.getidx(l_max, l, m)]
			sh_series_set(series, l, m, x)
			# reconstruct the missing negative-m coefficients
			if m:
				sh_series_set(series, l, -m, -x.conjugate() if (m & 1) else x.conjugate())
	return series


def healpy_alm_to_map(alms, l_max, m_max):
	"""
	Convert an a_{l,m} array and m_max integer to a pixel map.  Returns
	a 3-tuple containing three arrays (theta, phi, m).

	The combination

	>>> healpy_alm_to_map(*sh_series_to_healpy_alm(series))

	converts an sh_series object to a healpy pixel map.  The
	combination

	>>> healpy_alm_to_map(*read_fits_healpy_alm(filename))

	returns three arrays compatible with the return value of
	read_fits_healpy_map() but from a file containing a_{l,m}
	coefficients instead of the map itself.

	NOTE:  Raises ValueError if the (l_max, m_max) combination is not
	compatible with the requirements of the sphradiometer library (even
	if healpix allows it).
	"""
	if m_max not in (0, l_max):
		raise ValueError("unsupported (l_max, m_max) = (%d, %d)" % (l_max, m_max))

	m = healpy.alm2map(alms, int(2**math.ceil(math.log(l_max / 2., 2))), l_max, m_max)
	npix, = m.shape
	nside = healpy.npix2nside(npix)
	theta, phi = healpy.pixelfunc.pix2ang(nside, range(npix))
	return theta, phi, m


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


def read_fits_healpy_map(filename):
	m = healpy.fitsfunc.read_map(filename)
	npix, = m.shape
	nside = healpy.npix2nside(npix)
	theta, phi = healpy.pixelfunc.pix2ang(nside, range(npix))
	return theta, phi, m


def read_fits_healpy_alm(filename):
	"""
	Returns the tuple (alms, l_max, m_max), where alms is a numpy array
	containing the complex-valued a_{l,m} coefficients of the expansion
	of a real-valued function on the sphere, and l_max and m_max are
	the (integer) maximum l and m indexes for which there are
	coefficients.  The alms array is in healpix' native coefficient
	ordering.

	NOTE:  healpix does not support complex-valued functions on the
	sphere, so the negative m coefficients are not present.

	See also:  healpy_alm_to_map() to convert the tuple returned by
	this function into a pixel map.
	"""
	alms, m_max = healpy.fitsfunc.read_alm(filename, return_mmax = True)
	return alms, healpix_alm_l_max(alms, m_max), m_max
