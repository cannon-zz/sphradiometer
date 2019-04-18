# Copyright (C) 2006--2009,2019  Kipp C. Cannon
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


import re

from . import sphradiometer


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


def read_sh_series(f):
	pattern = re.compile(r"\((?P<l>\S+),(?P<m>\S+)\) = (?P<re>\S+) \+ I (?P<im>\S+)")
	lmax = 0
	polar = 1
	lines = []
	for line in f:
		l, m, r, i = re.search(pattern, line).groups()
		l, m, r, i = int(l), int(m), float(r), float(i)
		lines.append((l, m, r, i))
		lmax = max(lmax, l)
		polar &= m == 0
	if not lines:
		return radiometer.sh_series_new_zero(0, 0)
	series = radiometer.sh_series_new_zero(lmax, polar)
	for l, m, r, i in lines:
		radiometer.sh_series_set(series, l, m, complex(r, i))
	return series
