/*
 * Copyright (C) 2006,2009,2019  Kipp C. Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <sys/time.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/inject.h>
#include <sphradiometer/correlator.h>
#include <sphradiometer/diagnostics.h>


/*
 * ============================================================================
 *
 *                             Analytic Solutions
 *
 * ============================================================================
 */


/*
 * White Gaussian noise confined to Nyquist pass-band.
 */


static double sinc(double w)
{
	w *= M_PI;
	return w == 0.0 ? 1.0 : sin(w) / w;
}


struct __gaussian_noise_solution_data {
	double variance;
	double f0, delta_f;
	double baseline_d;
	double baseline_theta, baseline_phi;
	double cos_gamma_s;
};


static double __gaussian_noise_solution(double theta, double phi, void *Data)
{
	struct __gaussian_noise_solution_data *data = Data;

	/* cosine of the angle between the baseline and (theta, phi) */
	double cos_gamma = cos(theta) * cos(data->baseline_theta) + sin(theta) * sin(data->baseline_theta) * cos(phi - data->baseline_phi);

	double zeta = data->baseline_d * (cos_gamma - data->cos_gamma_s);

	/* formula is "spectral density * band width * ...", but the
	 * spectral density is the variance divided by the bandwidth, so we
	 * can just multiply by the variance */
	return data->variance * sinc(zeta * data->delta_f) * cos(2 * M_PI * zeta * data->f0);
}


static struct sh_series *gaussian_noise_solution(struct correlator_baseline *baseline, unsigned int l, double variance, double delta_t, double ra, double dec, double gmst)
{
	struct sh_series *series = sh_series_new(l, 0);
	double injection_theta = M_PI_2 - dec;
	double injection_phi = ra;
	struct __gaussian_noise_solution_data data = {
		/* variance */
		.variance = variance,
		/* pass band */
		.f0 = 0.25 / delta_t,
		.delta_f = 0.5 / delta_t,
		/* baseline description */
		.baseline_d = vector_magnitude(baseline->d),
		.baseline_theta = baseline->theta,
		.baseline_phi = baseline->phi + gmst,
		/* cosine of the angle between the baseline and the source */
		.cos_gamma_s = cos(injection_theta) * cos(baseline->theta) + sin(injection_theta) * sin(baseline->theta) * cos(injection_phi - (baseline->phi + gmst))
	};

	return sh_series_from_realfunc(series, __gaussian_noise_solution, &data);
}


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


static void usage(void)
{
	fprintf(stderr, "--sample-frequency	set the sample frequency\n");
	fprintf(stderr, "--trials	run this many trials\n");
}


struct options {
	int sample_frequency;
};


static struct options *command_line_options_new(void)
{
	struct options *options = malloc(sizeof(*options));

	/* defaults */
	*options = (struct options) {
		.sample_frequency = 512,
	};

	return options;
}


static void command_line_options_free(struct options *options)
{
	if(options) {
	}
	free(options);
}


struct options *parse_command_line(int argc, char *argv[])
{
	int c;
	int option_index;
	struct option long_options[] = {
		{"help",	0,	NULL,	'h'},
		{"sample-frequency",	required_argument,	NULL,	'f'},
		{NULL,	0,	NULL,	0}
	};
	struct options *options = command_line_options_new();

	if(!options)
		return NULL;

	do switch(c = getopt_long(argc, argv, "fht", long_options, &option_index)) {
	/* --sample-frequency */
	case 'f':
		options->sample_frequency = atoi(optarg);
		break;

	/* --help */
	case 'h':
		usage();
		exit(0);
		break;

	/* option sets a flag */
	case 0:
		break;

	/* end of arguments */
	case -1:
		break;

	/* unrecognized option */
	case '?':
		break;

	/* missing argument for an option */
	case ':':
		break;
	} while(c != -1);

	return options;
}


/*
 * ============================================================================
 *
 *                                   Tests
 *
 * ============================================================================
 */


static double dB(double x)
{
	return 10.0 * log10(fabs(x));
}


static int test1(gsl_rng *rng, struct instrument_array *instruments, double delta_t, double *gmst, int n_gmst)
{
	struct correlator_network_baselines *baselines = correlator_network_baselines_new(instruments);

	const int sky_l_max = correlator_network_l_max(baselines, delta_t);
	const int time_series_length = correlator_dump_interval(sky_l_max, 20) / delta_t;

	double injection[time_series_length];
	double *tseries[instrument_array_len(instruments)];
	complex double *fseries[instrument_array_len(instruments)];
	fftw_plan fftplans[instrument_array_len(instruments)];

	struct correlator_network_plan_fd *fdplans = correlator_network_plan_fd_new(baselines, time_series_length, delta_t);
	struct correlator_network_plan_td *tdplans = correlator_network_plan_td_new(baselines, delta_t);

	double *windows[baselines->n_baselines];

	struct sh_series *tdsky = sh_series_new_zero(sky_l_max, 0);
	struct sh_series *fdsky = sh_series_new_zero(sky_l_max, 0);

	const double variance = 1.0;

	int N_T, l_T, l_xi;

	int i, j;

	for(j = 0; j < instrument_array_len(instruments); j++) {
		tseries[j] = malloc(time_series_length * sizeof(**tseries));
		fseries[j] = malloc(time_series_length * sizeof(**fseries));
		fftplans[j] = correlator_tseries_to_fseries_plan(tseries[j], fseries[j], time_series_length);
		if(!tseries[j] || !fseries[j] || !fftplans[j]) {
			fprintf(stderr, "memory allocation failure\n");
			do {
				free(tseries[j]);
				free(fseries[j]);
				fftw_destroy_plan(fftplans[j]);
			} while(--j >= 0);
			return -1;
		}
	}
	for(j = 0; j < baselines->n_baselines; j++) {
		windows[j] = correlator_square_window_new(time_series_length - 2 * tdplans->plans[j]->transient, 0, 1.0 / (time_series_length - 2 * tdplans->plans[j]->transient));
		if(!windows[j]) {
			fprintf(stderr, "memory allocation failure\n");
			return -1;
		}
	}

	fprintf(stderr, "=== Parameters ===\n");
	fprintf(stderr, "=== Time-Domain Correlator ===\n");
	diagnostics_dump_correlator_plan_td_stats(stderr, tdplans->plans[0]);
	fprintf(stderr, "=== Frequency-Domain Correlator ===\n");
	diagnostics_dump_correlator_plan_fd_stats(stderr, fdplans->plans[0]);
	fprintf(stderr, "=== Other Numbers ===\n");
	fprintf(stderr, "l_{xi} = %d\n", sky_l_max);
	fprintf(stderr, "=== End Parameters ===\n");

	/* extract correlator parameters */

	N_T = tdplans->plans[0]->proj_a->n;
	l_T = tdplans->plans[0]->proj_a->l_max;
	l_xi = sky_l_max;

	/* loop over GMST */
	for(i = 0; i < n_gmst; i++) {
		struct sh_series *solution, *solution_lo;
		struct timeval tstart, tend;
		double runtime;

		fprintf(stderr, "=== Begin Test ===\nGMST = %g\n", gmst[i]);

		/* zero the instrument time series */
		for(j = 0; j < instrument_array_len(instruments); j++)
			memset(tseries[j], 0, time_series_length * sizeof(*tseries[j]));

		/* inject point source at vernal equinox */
		gaussian_white_noise(injection, time_series_length, variance, rng);
		for(j = 0; j < instrument_array_len(instruments); j++)
			inject_into(instrument_array_get(instruments, j), tseries[j], injection, time_series_length, delta_t, 0.0, 0.0, gmst[i]);

		/* compute integrated cross power */

		/* time domain */
		assert(correlator_network_integrate_power_td(tdsky, tseries, time_series_length, windows, tdplans) != NULL);
		sh_series_rotate_z(tdsky, tdsky, gmst[i]);

		/* record the start time */
		gettimeofday(&tstart, NULL);

		/* frequency domain */
		for(j = 0; j < instrument_array_len(instruments); j++)
			correlator_tseries_to_fseries(tseries[j], fseries[j], time_series_length, fftplans[j]);
		assert(correlator_network_integrate_power_fd(fdsky, fseries, fdplans) != NULL);
		sh_series_rotate_z(fdsky, fdsky, gmst[i]);

		/* record the end time */
		gettimeofday(&tend, NULL);

		/* compute the execution speed */
		runtime = (tend.tv_sec - tstart.tv_sec) + 1e-6 * (tend.tv_usec - tstart.tv_usec);
		fprintf(stderr, "=== Begin Speed ===\n");
		fprintf(stderr, "%d samples analyzed by FD correlator in %g seconds = %.3g samples/second\n", time_series_length, runtime, time_series_length / runtime);
		fprintf(stderr, "=== End Speed ===\n");

		/* compute exact-ish solution for source at vernal equinox.
		 * we still compute it to some finite harmonic order, but
		 * much higher than what the correlator will compute, so
		 * the residual errors in this solution should not
		 * contribute significantly to the overall error budget */
		solution = gaussian_noise_solution(baselines->baselines[0], 2 * sky_l_max, variance, delta_t, 0.0, 0.0, gmst[i]);
		for(j = 1; j < baselines->n_baselines; j++) {
			struct sh_series *tmp = gaussian_noise_solution(baselines->baselines[j], 2 * sky_l_max, variance, delta_t, 0.0, 0.0, gmst[i]);
			sh_series_add(solution, 1.0, tmp);
			sh_series_free(tmp);
		}
		sh_series_scale(solution, 1. / baselines->n_baselines);

		/* compute a version of the exact solution truncated to the
		 * harmonic order of the correlator.  by comparing this to
		 * the other, more exact, solution, we measure how much
		 * error is coming from the truncation of the sky to finite
		 * order.  any additional error observed in the correlator
		 * output is from the correlator. */
		solution_lo = sh_series_resize(sh_series_copy(solution), sky_l_max);

		/* solution output */
		{
		char filename[100];
		sprintf(filename, "tests_tdaverage_%d_%d_%d_%020.17f.fits", N_T, l_T, l_xi, gmst[i]);
		assert(sh_series_write_healpix_map(tdsky, filename) == 0);
		sprintf(filename, "tests_fdaverage_%d_%d_%d_%020.17f.fits", N_T, l_T, l_xi, gmst[i]);
		assert(sh_series_write_healpix_map(fdsky, filename) == 0);
		sprintf(filename, "tests_exact_%d_%d_%d_%020.17f.fits", N_T, l_T, l_xi, gmst[i]);
		assert(sh_series_write_healpix_map(solution, filename) == 0);
		}

		/* compute normalized RMS error */
		fprintf(stderr, "=== Begin Errors ===\n");
		fprintf(stderr, "RMS (exact - order-reduced exact) = %g dB\n", dB(diagnostics_rms_error(solution_lo, solution) / variance));
		fprintf(stderr, "=== Time-Domain Correlator Errors ===\n");
		fprintf(stderr, "RMS (TD - exact) = %g dB\n", dB(diagnostics_rms_error(tdsky, solution) / variance));
		fprintf(stderr, "RMS (TD - order-reduced exact) = %g dB\n", dB(diagnostics_rms_error(tdsky, solution_lo) / variance));
		fprintf(stderr, "TD fractional error @ source = %g dB\n", dB((creal(sh_series_eval(solution, 0.0, 0.0)) - creal(sh_series_eval(tdsky, 0.0, 0.0))) / variance));

		fprintf(stderr, "=== Frequency-Domain Correlator Errors ===\n");
		fprintf(stderr, "fractional RMS (FD - exact) = %g dB\n", dB(diagnostics_rms_error(fdsky, solution) / variance));
		fprintf(stderr, "fractional RMS (FD - order-reduced exact) = %g dB\n", dB(diagnostics_rms_error(fdsky, solution_lo) / variance));
		fprintf(stderr, "FD fractional error @ source = %g dB\n", dB((creal(sh_series_eval(solution, 0.0, 0.0)) - creal(sh_series_eval(fdsky, 0.0, 0.0))) / variance));
		fprintf(stderr, "fractional RMS (TD - FD) = %g dB\n", dB(diagnostics_rms_error(fdsky, tdsky) / variance));
		fprintf(stderr, "=== End Errors ===\n");

		if(dB(diagnostics_rms_error(tdsky, solution_lo) / variance) > -21. || dB(diagnostics_rms_error(fdsky, solution_lo) / variance) > -21.) {
			fprintf(stderr, "accuracy failure.\n");
			return -1;
		}

		/* prepare for next test */
		sh_series_free(solution);
		sh_series_free(solution_lo);

		fprintf(stderr, "=== End Test ===\n");
	}

	/*
	 * Clean up
	 */

	for(j = 0; j < instrument_array_len(instruments); j++) {
		free(tseries[j]);
		free(fseries[j]);
		fftw_destroy_plan(fftplans[j]);
	}
	for(j = 0; j < baselines->n_baselines; j++)
		free(windows[j]);
	correlator_network_plan_fd_free(fdplans);
	correlator_network_plan_td_free(tdplans);
	correlator_network_baselines_free(baselines);
	sh_series_free(tdsky);
	sh_series_free(fdsky);

	return 0;
}


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	struct options *options = parse_command_line(argc, argv);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);

	struct instrument_array *instruments = instrument_array_new(0);
	instrument_array_append(instruments, instrument_new_from_name("H1"));
	instrument_array_append(instruments, instrument_new_from_name("L1"));
	instrument_array_append(instruments, instrument_new_from_name("V1"));
	instrument_array_append(instruments, instrument_new_from_name("K1"));

	const double delta_t = 1.0 / options->sample_frequency;

	const int n_gmst = 3;
	double gmst[n_gmst];
	int i;

	gsl_rng_set(rng, time(NULL));

	for(i = 0; i < n_gmst; i++)
		gmst[i] = i * (M_PI / 4);

	assert(test1(rng, instruments, delta_t, gmst, n_gmst) == 0);


	/*
	 * Clean up
	 */

	instrument_array_free(instruments);
	command_line_options_free(options);
	gsl_rng_free(rng);

	return 0;
}
