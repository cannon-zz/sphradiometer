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
	if m_max == 0:
		l_max, = alms.shape
	elif alms.shape != (((m_max + 1) * (m_max + 2)) // 2,):
		# for mmax == lmax, when no negative m coefficients are
		# present, the total number of coefficients must be related
		# to m_max as above.  if they aren't, then this is a kind
		# of healpix alm set that we don't support
		raise ValueError("unsupported size")
	else:
		l_max = m_max
	return alms, l_max, m_max


def healpy_alm_to_map(alms, l_max, m_max):
	"""
	Convert an a_{l,m} array and m_max integer to a pixel map.  The
	combination

	>>> healpy_alm_to_map(*read_fits_healpy_alm(filename))

	returns three arrays compatible with the return value of
	read_fits_healpy_map() but from a file containing a_{l,m}
	coefficients instead of the map itself.
	"""
	if m_max not in (0, l_max):
		raise ValueError("bad m_max")

	m = healpy.alm2map(alms, int(2**math.ceil(math.log(l_max / 2., 2))), l_max, m_max)
	npix, = m.shape
	nside = healpy.npix2nside(npix)
	theta, phi = healpy.pixelfunc.pix2ang(nside, range(npix))
	return theta, phi, m
