# Copyright (C) 2006--2009  Kipp C. Cannon
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
import matplotlib
matplotlib.rcParams.update({
	"text.usetex": True,
	"figure.dpi" : 300,
	"savefig.dpi" : 300
})
from matplotlib import cm
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
import numpy

from . import sphradiometer


#
# =============================================================================
#
#                                  Some Math
#
# =============================================================================
#


def rad_to_hms(a):
	# wrap a into [0, 2 pi)
	a %= 2 * math.pi
	# convert to seconds
	a *= 12 * 60 * 60 / math.pi
	# split and return (h, m, s)
	h = int(a / 3600)
	a %= 3600
	return h, int(a / 60), a % 60


def rad_to_hms_str(a):
	return "%02d h %02d m %.2f s" % rad_to_hms(a)


def rad_to_dms(a):
	# wrap a into [-pi, pi)
	a = (a + math.pi) % (2 * math.pi) - math.pi
	# convert to arcseconds
	a *= 180 * 60 * 60 / math.pi
	# split and return (d, m, s)
	d = int(a / 3600)
	a = abs(a) % 3600
	return d, int(a / 60), a % 60


def rad_to_dms_str(a):
	return "%+02do %02d' %.2f\"" % rad_to_dms(a)


#
# =============================================================================
#
#                                    Plots
#
# =============================================================================
#


class XYSlicePlot(object):
	def __init__(self, axes, pixels_per_pi = 4096):
		self.axes = axes

		# sample co-ordinates uniformly distributed around equator
		self.phi = numpy.linspace(0., 2. * math.pi, 2 * pixels_per_pi + 1)
		self.theta = numpy.zeros((len(self.phi),), dtype = "double") + math.pi / 2.

		# put a dummy point on the axes to force a scale
		self.axes.plot((0,), (2.2,))

	def yvals_from_sh_series(self, series):
		yvals = numpy.zeros((len(self.phi),), dtype = "double")
		for i in range(len(self.phi)):
			yvals[i] = sphradiometer.sh_series_eval(series, self.theta[i], self.phi[i]).real
		return yvals

	def yvals_from_healpix_map(self, m):
		yvals = numpy.zeros((len(self.phi),), dtype = "double")
		for i in range(len(self.phi)):
			yvals[i] = healpy.pixelfunc.get_interp_val(m,  self.theta[i], self.phi[i])
		return yvals

	def plot_yvals(self, yvals, *args, **kwargs):
		self.axes.plot(self.phi, yvals, *args, **kwargs)


class SkyPlot(object):
	def __init__(self, axes, pixels_per_pi = 256):
		self.map = Basemap(projection="moll", rsphere = 1, lon_0 = 0, ax = axes)
		self.map.drawparallels(numpy.linspace(-90., +90, 13), linewidth = 0.3)
		self.map.drawmeridians(numpy.linspace(-180., +180., 25), linewidth = 0.3)

		self.ra = numpy.linspace(-math.pi, math.pi, 2 * pixels_per_pi + 1)
		self.dec = numpy.linspace(-math.pi / 2., math.pi / 2., pixels_per_pi + 1)

	def sh_series_contourf(self, series, cmap = cm.gray_r, **kwargs):
		samples = numpy.zeros((len(self.ra), len(self.dec)), dtype = "double")
		# use -ra for the azimuthal co-ordinate when evaluating the
		# function, and +ra for the values passed to meshgrid() so
		# that the sky is flipped around the right way.  flipping
		# the sky by, instead, passing -ra to meshgrid() produces a
		# matplotlib warning about a non-monotonically increasing x
		# co-ordinate
		for i, ra in enumerate(self.ra):
			for j, dec in enumerate(self.dec):
				samples[i, j] = sphradiometer.sh_series_eval(series, math.pi / 2. - dec, -ra).real
		x, y = self.map(*numpy.meshgrid(self.ra * 180. / math.pi, self.dec * 180. / math.pi))
		self.map.contourf(x, y, numpy.transpose(samples), cmap = cmap, **kwargs)

	def healpix_map_contourf(self, m, cmap = cm.gray_r, **kwargs):
		samples = numpy.zeros((len(self.ra), len(self.dec)), dtype = "double")
		# use -ra for the azimuthal co-ordinate when evaluating the
		# function, and +ra for the values passed to meshgrid() so
		# that the sky is flipped around the right way.  flipping
		# the sky by, instead, passing -ra to meshgrid() produces a
		# matplotlib warning about a non-monotonically increasing x
		# co-ordinate
		for i, ra in enumerate(self.ra):
			for j, dec in enumerate(self.dec):
				samples[i, j] = healpy.pixelfunc.get_interp_val(m,  math.pi / 2. - dec, -ra)
		x, y = self.map(*numpy.meshgrid(self.ra * 180. / math.pi, self.dec * 180. / math.pi))
		self.map.contourf(x, y, numpy.transpose(samples), cmap = cmap, **kwargs)

