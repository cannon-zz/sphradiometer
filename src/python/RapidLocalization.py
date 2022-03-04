import numpy as np
from os.path import exists

from . import sphradiometer as sph


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
#                                 Localization
#
# ==============================================================================
#


class RapidLocalization(object):
	def __init__(self, instruments, psds, length, deltaT):
		"""
		Parameter
		---------
		instruments : list
			list of observing detectors.  e.g., ["H1", "L1", "V1"]
		psds : <Swig Object of type 'double **'>
		length : int
			legth of SNR time series used for localization.
		deltaT : float
			bin width of SNR time series used for localization.

		NOTE
		----
		If the length in reading-process of the precalculated objects
		is NOT equal to that in writing-process, then the precalculated
		objects have to be generated again.
		"""
		# preamble
		self.inst_array = sph.instrument_array_new(0)
		for ifo in instruments:
			sph.instrument_array_append(self.inst_array, sph.instrument_new_from_name(ifo))

		# prepare precalculated objects
		self.baselines = sph.correlator_network_baselines_new(self.inst_array)
		if len(instruments) > 1:
			# multi detector case
			self.logprior = sph.sh_series_log_uniformsky_prior(sph.correlator_network_l_max(self.baselines, deltaT))

			self.fdplansp = sph.correlator_network_plan_fd_new(self.baselines, length, deltaT)
			self.fdplansn = sph.correlator_network_plan_fd_copy(self.fdplansp)
			sph.correlator_network_plan_mult_by_projection(self.fdplansp, +1, 0, psds)
			sph.correlator_network_plan_mult_by_projection(self.fdplansn, -1, 0, psds)

			self.fdautoplanp = sph.autocorrelator_network_plan_fd_new(self.inst_array, +1, 0, psds, int(sph.pick_ith_correlator_plan_fd(self.fdplansp.plans, 0).delay_product.n), self.logprior.l_max)
			self.fdautoplann = sph.autocorrelator_network_plan_fd_new(self.inst_array, -1, 0, psds, int(sph.pick_ith_correlator_plan_fd(self.fdplansp.plans, 0).delay_product.n), self.logprior.l_max)

		else:
			# single detector case
			self.logprior = sph.sh_series_log_uniformsky_prior(sph.correlator_baseline_power_l_max_naive(length * deltaT, deltaT))


	def sphcoeff(self, skyp, skyn, seriesp, nseriesp):
		"""
		Parameter
		---------
		skyp : a coeff_series instance
			object to store a localization result for \\beta = +1
		skyn : a coeff_series instance
			object to store a localization result for \\beta = -1
		seriesp : <Swig Object of type 'COMPLEX16TimeSeries **'>
			stuff of observed SNR time series for all detectors
		seriesn : <Swig Object of type 'COMPLEX16TimeSeries **'>
			stuff of template auto-correlations on time domain for
			all detectors

		NOTE
		----
		seriesp & seriesn have to be pre-processed by
		preprocess_SNRTimeSeries()
		"""
		if sph.instrument_array_len(self.inst_array) > 1:
			# multi detector case
			sph.generate_alm_skys(skyp.coeff, skyn.coeff, self.fdplansp, self.fdplansn, self.fdautoplanp, self.fdautoplann, seriesp, nseriesp, self.logprior)
		else:
			# single detector case
			sph.sh_seriesp_assign(skyp.coeff, self.logprior)
			sph.sh_seriesp_assign(skyn.coeff, self.logprior)

	def write(self, precalc_path):
		"""
		Parameter
		---------
		precalc_path : str
			path (name) of used precalculated objects.
		"""
		sph.make_precalc_directories(precalc_path, sph.instrument_array_len(self.inst_array))
		if sph.instrument_array_len(self.inst_array) > 1:
			# (n > 1) detector case
			sph.write_precalc_logprior(self.logprior, precalc_path)
			sph.write_precalc_correlator_network_plan_fd(self.fdplansp, self.fdplansn, precalc_path)
			sph.write_precalc_autocorrelator_network_plan_fd(self.fdautoplanp, self.fdautoplann, precalc_path)
		else:
			# 1 detector case
			sph.write_precalc_logprior(self.logprior, precalc_path)

	@classmethod
	def read(cls, precalc_path, instruments, length):
		"""
		Parameter
		---------
		precalc_path : str
			path (name) of used precalculated objects.
		"""
		# FIXME:  it is a design error that the data cannot be
		# loaded off disk unless the calling code knows the
		# instrument list, SNR time series length, etc..  imagine a
		# word processor being unable to load a document unless you
		# can remember the title of chapter 2.  what is written to
		# disk must be self-describing, so that a program can
		# simply load the data.

		self = cls.__new__(cls)

		# preamble
		self.inst_array = sph.instrument_array_new(0)
		for ifo in instruments:
			sph.instrument_array_append(self.inst_array, sph.instrument_new_from_name(ifo))
		self.baselines = sph.correlator_network_baselines_new(self.inst_array)

		if sph.instrument_array_len(self.inst_array) > 1:
			# (n > 1) detector case
			self.logprior = sph.read_precalc_logprior(precalc_path)

			self.fdplansp = sph.new_correlator_network_plan_fdp()
			self.fdplansn = sph.new_correlator_network_plan_fdp()
			sph.read_precalc_correlator_network_plan_fd(self.fdplansp, self.fdplansn, self.inst_array, length, precalc_path)

			self.fdautoplanp = sph.new_autocorrelator_network_plan_fdp()
			self.fdautoplann = sph.new_autocorrelator_network_plan_fdp()
			sph.read_precalc_autocorrelator_network_plan_fd(self.fdautoplanp, self.fdautoplann, self.inst_array, length, precalc_path)
		else:
			# 1 detector case
			self.logprior = sph.read_precalc_logprior(precalc_path)

		return self

	def __del__(self):
		sph.correlator_network_baselines_free(self.baselines)
		sph.sh_series_free(self.logprior)
		if sph.instrument_array_len(self.inst_array) > 1:
			sph.correlator_network_plan_fd_free(self.fdplansp)
			sph.correlator_network_plan_fd_free(self.fdplansn)
			sph.autocorrelator_network_plan_fd_free(self.fdautoplanp)
			sph.autocorrelator_network_plan_fd_free(self.fdautoplann)
		sph.instrument_array_free(self.inst_array)
