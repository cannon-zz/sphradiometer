from glob import glob
from itertools import combinations
import numpy as np
import os

from lal import series as lalseries
from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import param as ligolw_param
from ligo.lw import utils as ligolw_utils

from . import sphradiometer as sph


@lsctables.use_in
class ContentHandler(lalseries.PSDContentHandler):
	pass


#
# ==============================================================================
#
#                           Masked Memory Operation
#
# ==============================================================================
#


class coeff_series(object):
	def __init__(self, ndet):
		"""
		Parameter
		---------
		ndet : int
			the number of observing detectors
		"""
		self.ndet = ndet
		if ndet > 1:
			self.coeff = sph.new_sh_seriespp()
		else:
			self.coeff = sph.new_sh_seriesp()

	def get(self):
		if self.ndet > 1:
			return sph.sh_seriespp_value(self.coeff)
		else:
			return self.coeff

	def __del__(self):
		if self.ndet > 1:
			sph.sh_series_free(sph.sh_seriespp_value(self.coeff))
			sph.delete_sh_seriespp(self.coeff)
		else:
			sph.delete_sh_seriesp(self.coeff)


#
# ==============================================================================
#
#                                    PSDs
#
# ==============================================================================
#


class flat_psd(object):
	def __init__(self, length):
		"""
		Parameter
		---------
		length : int
			psd length
		"""
		self.psd = sph.new_double_array(length)
		for i in range(length):
			sph.double_array_setitem(self.psd, i, 1)

	def __del__(self):
		sph.delete_double_array(self.psd)


def transpose_psd_array(psds, length):
	"""
	Parameter
	---------
	psds : dict of e.g. flat_psd instances
		key is detector prefix
	length : int
		psd length
	"""
	instruments = sorted(psds)

	psds_array = sph.new_doublep_array(len(psds))
	for i, ifo in enumerate(instruments):
		sph.doublep_array_setitem(psds_array, i, psds[ifo].psd)
	tpsds_array = sph.transpose_matrix(psds_array, len(psds), length)
	sph.delete_doublep_array(psds_array)

	return tpsds_array


#
# ==============================================================================
#
#                                   Time Series
#
# ==============================================================================
#


class create_sph_COMPLEX16TimeSeries_wrap(object):
	def __init__(self, tseries):
		"""
		Parameter
		---------
		tseries : <Swig Object of type 'COMPLEX8TimeSeries *'>
		"""
		data = tseries.data.data
		data_ = sph.new_double_complex_array(len(data))
		for i, x in enumerate(data):
			sph.double_complex_array_setitem(data_, i, complex(x))

		self.series = sph.new_COMPLEX16TimeSeriesp()
		sph.create_sph_COMPLEX16TimeSeries(
			self.series,
			tseries.name,
			tseries.epoch.gpsSeconds,
			tseries.epoch.gpsNanoSeconds,
			tseries.f0,
			tseries.deltaT,
			None,
			len(data),
			data_
		)

		sph.delete_double_complex_array(data_)

	def __del__(self):
		sph.free_SNRTimeSeries(self.series)


class convert_TimeSeries2Sequence_wrap(object):
	def __init__(self, series):
		"""
		Parameter
		---------
		series : create_sph_COMPLEX16TimeSeries_wrap instances
		"""
		self.series = sph.convert_TimeSeries2Sequence(series.series)

	def __del__(self):
		sph.free_SNRSequence(self.series)


#
# ==============================================================================
#
#                                 Localization
#
# ==============================================================================
#


class RapidLocalization_(object):
	def __init__(self, psds, precalc_length, deltaT):
		"""
		Parameter
		---------
		psds : dict of e.g. flat_psd instances
			psds for all detectors
		precalc_length : int
			legth of SNR time series used for localization.
		deltaT : float
			bin width of SNR time series used for localization.

		NOTE
		----
		If the precalc_length in reading-process of the precalculated
		objects is NOT equal to that in writing-process, then the
		precalculated objects have to be generated again.
		"""
		# preamble
		self.precalc_length = precalc_length
		self.instruments = sorted(psds)
		self.inst_array = sph.instrument_array_new(0)
		for ifo in self.instruments:
			sph.instrument_array_append(self.inst_array, sph.instrument_new_from_name(ifo))
		psd_array = transpose_psd_array(psds, self.precalc_length)

		# prepare precalculated objects
		self.baselines = sph.correlator_network_baselines_new(self.inst_array)
		if len(self.instruments) > 1:
			# multi detector case
			self.logprior = sph.sh_series_log_uniformsky_prior(sph.correlator_network_l_max(self.baselines, deltaT))

			self.fdplansp = sph.correlator_network_plan_fd_new(self.baselines, self.precalc_length, deltaT)
			self.fdplansn = sph.correlator_network_plan_fd_copy(self.fdplansp)
			sph.correlator_network_plan_mult_by_projection(self.fdplansp, +1, 0, psd_array)
			sph.correlator_network_plan_mult_by_projection(self.fdplansn, -1, 0, psd_array)

			self.fdautoplanp = sph.autocorrelator_network_plan_fd_new(self.inst_array, +1, 0, psd_array, int(sph.pick_ith_correlator_plan_fd(self.fdplansp.plans, 0).delay_product.n), self.logprior.l_max)
			self.fdautoplann = sph.autocorrelator_network_plan_fd_new(self.inst_array, -1, 0, psd_array, int(sph.pick_ith_correlator_plan_fd(self.fdplansp.plans, 0).delay_product.n), self.logprior.l_max)

		else:
			# single detector case
			self.logprior = sph.sh_series_log_uniformsky_prior(sph.correlator_baseline_power_l_max_naive(self.precalc_length * deltaT, deltaT))

		# free
		for i in range(self.precalc_length):
			sph.delete_double_array(sph.doublep_array_getitem(psd_array, i))
		sph.delete_doublep_array(psd_array)


	def sphcoeff(self, snr, aut):
		"""
		Parameter
		---------
		snr : dict of <Swig Object of type 'COMPLEX16TimeSeries *'>
			stuff of observed SNR time series for all detectors
		aut : dict of <Swig Object of type 'COMPLEX16TimeSeries *'>
			stuff of template auto-correlations on time domain for
			all detectors
		Returns
		-------
		output : tuple of coeff_series instances
			object to store a localization result for \\beta = +/-1
		"""
		# construct **COMPLEX16_array
		seriesp = sph.new_COMPLEX16TimeSeries_array(len(snr))
		nseriesp = sph.new_COMPLEX16Sequence_array(len(aut))
		for i, ifo in enumerate(self.instruments):
			sph.COMPLEX16TimeSeries_array_setitem(seriesp, i, snr[ifo].series)
			sph.COMPLEX16Sequence_array_setitem(nseriesp, i, aut[ifo].series)
		sph.preprocess_SNRTimeSeries(seriesp, nseriesp, len(self.instruments))

		# consistency check
		assert self.instruments == sorted(snr.keys()), "instruments are inconsistent"
		if len(snr) > 1:
			assert self.precalc_length == sph.pick_length_from_COMPLEX16TimeSeries(sph.COMPLEX16TimeSeries_array_getitem(seriesp, 0)), "SNR time series length is inconsistent"

		# store sph coefficient series result
		skyp = coeff_series(len(self.instruments))
		skyn = coeff_series(len(self.instruments))
		if sph.instrument_array_len(self.inst_array) > 1:
			# multi detector case
			sph.generate_alm_skys(skyp.coeff, skyn.coeff, self.fdplansp, self.fdplansn, self.fdautoplanp, self.fdautoplann, seriesp, nseriesp, self.logprior)
		else:
			# single detector case
			sph.sh_seriesp_assign(skyp.coeff, self.logprior)
			sph.sh_seriesp_assign(skyn.coeff, self.logprior)

		# free **COMPLEX16_array
		sph.delete_COMPLEX16TimeSeries_array(seriesp)
		sph.delete_COMPLEX16Sequence_array(nseriesp)

		return skyp, skyn

	def write(self, precalc_path):
		"""
		Parameter
		---------
		precalc_path : str
			path (name) of used precalculated objects.
		"""
		sph.make_precalc_directories(precalc_path, sph.instrument_array_len(self.inst_array))
		np.save(os.path.join(precalc_path, "instruments.npy"), self.instruments)
		sph.write_precalc_time_series_length(self.precalc_length, precalc_path)
		if sph.instrument_array_len(self.inst_array) > 1:
			# (n > 1) detector case
			sph.write_precalc_logprior(self.logprior, precalc_path)
			sph.write_precalc_correlator_network_plan_fd(self.fdplansp, self.fdplansn, precalc_path)
			sph.write_precalc_autocorrelator_network_plan_fd(self.fdautoplanp, self.fdautoplann, precalc_path)
		else:
			# 1 detector case
			sph.write_precalc_logprior(self.logprior, precalc_path)

	@classmethod
	def read(cls, precalc_path):
		"""
		Parameter
		---------
		precalc_path : str
			path (name) of used precalculated objects.
		"""
		self = cls.__new__(cls)

		# preamble
		self.instruments = list(np.load(os.path.join(precalc_path, "instruments.npy")))
		self.inst_array = sph.instrument_array_new(0)
		for ifo in self.instruments:
			sph.instrument_array_append(self.inst_array, sph.instrument_new_from_name(ifo))
		self.baselines = sph.correlator_network_baselines_new(self.inst_array)

		self.logprior = sph.read_precalc_logprior(precalc_path)
		if sph.instrument_array_len(self.inst_array) > 1:
			# (n > 1) detector case
			length = sph.new_uintp()
			sph.read_precalc_time_series_length(length, precalc_path)
			self.precalc_length = sph.uintp_value(length)
			sph.delete_uintp(length)

			self.fdplansp = sph.new_correlator_network_plan_fdp()
			self.fdplansn = sph.new_correlator_network_plan_fdp()
			sph.read_precalc_correlator_network_plan_fd(self.fdplansp, self.fdplansn, self.inst_array, self.precalc_length, precalc_path)

			self.fdautoplanp = sph.new_autocorrelator_network_plan_fdp()
			self.fdautoplann = sph.new_autocorrelator_network_plan_fdp()
			sph.read_precalc_autocorrelator_network_plan_fd(self.fdautoplanp, self.fdautoplann, self.inst_array, self.precalc_length, precalc_path)

		return self

	def __del__(self):
		sph.sh_series_free(self.logprior)
		if sph.instrument_array_len(self.inst_array) > 1:
			sph.correlator_network_plan_fd_free(self.fdplansp)
			sph.correlator_network_plan_fd_free(self.fdplansn)
			sph.autocorrelator_network_plan_fd_free(self.fdautoplanp)
			sph.autocorrelator_network_plan_fd_free(self.fdautoplann)
		sph.instrument_array_free(self.inst_array)
		sph.correlator_network_baselines_free(self.baselines)


class RapidLocalization(object):
	def __init__(self, psds, precalc_length, deltaT):
		"""
		Parameter
		---------
		psds : dict of e.g. flat_psd instances
			psds for all detectors
		precalc_length : int
			legth of SNR time series used for localization.
		deltaT : float
			bin width of SNR time series used for localization.

		NOTE
		----
		If the precalc_length in reading-process of the precalculated
		objects is NOT equal to that in writing-process, then the
		precalculated objects have to be generated again.
		"""
		self.instruments = sorted(psds)
		self.precalcs = {}
		for n in range(1, len(self.instruments) + 1):
			for c in combinations(self.instruments, n):
				print("initialize for", "".join(c))
				psds_ = {ifo : psds[ifo] for ifo in c}
				self.precalcs["".join(c)] = RapidLocalization_(
					psds_,
					precalc_length,
					deltaT
				)

	def sphcoeff(self, snr, aut):
		"""
		Parameter
		---------
		snr : dict of <Swig Object of type 'COMPLEX16TimeSeries *'>
			stuff of observed SNR time series for all detectors
		aut : dict of <Swig Object of type 'COMPLEX16TimeSeries *'>
			stuff of template auto-correlations on time domain for
			all detectors
		Returns
		-------
		output : tuple of coeff_series instances
			object to store a localization result for \\beta = +/-1
		"""
		return self.precalcs["".join(sorted(snr))].sphcoeff(snr, aut)

	def write(self, precalc_path):
		"""
		Parameter
		---------
		precalc_path : str
			path (name) of used precalculated objects.
		"""
		os.mkdir(precalc_path)
		for n in range(1, len(self.instruments) + 1):
			for c in combinations(self.instruments, n):
				name = "".join(c)
				print("write for", name)
				self.precalcs[name].write(os.path.join(precalc_path, name))

	@classmethod
	def read(cls, precalc_path):
		"""
		Parameter
		---------
		precalc_path : str
			path (name) of used precalculated objects.
		"""
		self = cls.__new__(cls)

		self.precalcs = {}
		for prei in [os.path.basename(d) for d in glob(os.path.join(precalc_path, "*"))]:
			print("read for", prei)
			self.precalcs[prei] = RapidLocalization_.read(os.path.join(precalc_path, prei))
		return self


#
# ==============================================================================
#
#                                    Main
#
# ==============================================================================
#


if __name__ == "__main__":
	#
	# prepare data
	#

	print("preamble")
	coinc_file = "coinc.xml"
	precalc_path = "precalcs"
	print("assigned precalculated object is", precalc_path)

	snr = {}
	aut = {}
	xmldoc = ligolw_utils.load_filename(coinc_file, contenthandler=ContentHandler, verbose=True)
	sngl_inspiral_index = dict((row.event_id, row) for row in lsctables.SnglInspiralTable.get_table(xmldoc))
	for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName):
		if elem.hasAttribute("Name") and elem.Name == "COMPLEX8TimeSeries":
			sngl_inspiral = sngl_inspiral_index[int(ligolw_param.get_pyvalue(elem, "event_id"))]
			snr[sngl_inspiral.ifo] = create_sph_COMPLEX16TimeSeries_wrap(lalseries.parse_COMPLEX8TimeSeries(elem))
			# Currently nseriesp is dummy information, so that seriesp is set
			aut[sngl_inspiral.ifo] = convert_TimeSeries2Sequence_wrap(snr[sngl_inspiral.ifo])
	instruments = sorted(snr.keys())

	# Currently psd is dummy information, so that it's set one
	precalc_length = sph.precalculated_TimeSeries_length(
		sph.pick_length_from_COMPLEX16TimeSeries(snr[instruments[0]].series),
		sph.pick_deltaT_from_COMPLEX16TimeSeries(snr[instruments[0]].series)
	)
	psds = {}
	for ifo in instruments:
		psds[ifo] = flat_psd(precalc_length)


	#
	# Localize
	#


	print("prepare precalculated objects")
	if os.path.exists(precalc_path):
		print("read objects")
		rapidloc = RapidLocalization.read(precalc_path)
	else:
		print("make objects")
		rapidloc = RapidLocalization(
			psds,
			precalc_length,
			sph.pick_deltaT_from_COMPLEX16TimeSeries(snr[instruments[0]].series)
		)
		print("write objects")
		rapidloc.write(precalc_path)

	# sph coeff series for \beta = +1 & -1
	print("calc")
	skyp, skyn = rapidloc.sphcoeff(snr, aut)

	# save coeff series
	print("save")
	sph.sh_series_write_healpix_alm(skyp.get(), "coeffp.fits")
	sph.sh_series_write_healpix_alm(skyn.get(), "coeffn.fits")


	#
	# free
	#


	print("free")
	del snr
	del aut
	del psds
	del rapidloc
	del skyp
	del skyn
