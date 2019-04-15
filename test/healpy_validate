#!/usr/bin/env python3

import sys
from sphradiometer.sphradiometer import sh_series_new, sh_series_lmoffset
from healpy.sphtfunc import Alm

series = sh_series_new(64, 0)
for l in range(series.l_max + 1):
	# FIXME:  m < 0 is messed up
	for m in range(0, l + 1):
		if sh_series_lmoffset(series, l, m) != Alm.getidx(series.l_max, l, m):
			print("for (l, m) = (%d, %d): sphradiometer's index s %d, healpy's index is %d" % (l, m, sh_series_lmoffset(series, l, m), Alm.getidx(series.l_max, l, m)))
			sys.exit(1)
