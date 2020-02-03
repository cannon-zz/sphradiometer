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
	lines = []
	for n, line in enumerate(f, 1):
		# remove leading and trailing whitespace
		line = line.strip()
		# skip comments
		if line[0] == "#":
			continue
		# parse the line
		match = re.search(pattern, line)
		if match is None:
			raise ValueError("line %d: unrecognized format" % n)
		l, m, r, i = match.groups()
		try:
			l = int(l)
			m = int(m)
			if l < 0 or abs(m) > l:
				raise ValueError("invalid l, m: %d, %d" % (l, m))
			lines.append((l, m, complex(float(r), float(i))))
		except ValueError as e:
			raise ValueError("line %d: %s" % (n, str(e)))
	if not lines:
		return sphradiometer.sh_series_new_zero(0, 0)
	lmax = max(line[0] for line in lines)
	polar = not any(line[1] for line in lines)
	series = sphradiometer.sh_series_new_zero(lmax, polar)
	for l, m, a in lines:
		sphradiometer.sh_series_set(series, l, m, a)
	return series
