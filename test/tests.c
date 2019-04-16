/*
 * Copyright (C) 2006  Kipp C. Cannon
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
	double injection_phi = -ra;
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
	int trials;
};


static struct options *command_line_options_new(void)
{
	struct options *options = malloc(sizeof(*options));

	/* defaults */
	*options = (struct options) {
		.sample_frequency = 2048,
		.trials = 5,
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
		{"trials",	required_argument,	NULL,	't'},
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

	/* --trials */
	case 't':
		options->trials = atoi(optarg);
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
 *                              sh_series Tests
 *
 * ============================================================================
 */


/*
 * Test projection.  The Y_{lm} test functions have been copied from
 * Wikipedia.
 */


static complex double Y_00_00(double theta, double phi, void *nul)
{
	return 1.0 / 2 * sqrt(1 / M_PI);
}


static complex double Y_07_01(double theta, double phi, void *nul)
{
	return -1.0 / 64 * sqrt(105 / (2 * M_PI)) * cexp(I * phi) * sin(theta) * (429 * pow(cos(theta), 6) - 495 * pow(cos(theta), 4) + 135 * pow(cos(theta), 2) - 5);
}


static complex double Y_08_n7(double theta, double phi, void *nul)
{
	return 3.0 / 64 * sqrt(12155 / (2 * M_PI)) * cexp(I * -7 * phi) * pow(sin(theta), 7) * cos(theta);
}


static complex double Y_10_03(double theta, double phi, void *nul)
{
	return -3.0 / 256 * sqrt(5005 / M_PI) * cexp(I * 3 * phi) * pow(sin(theta), 3) * (323 * pow(cos(theta), 7) - 357 * pow(cos(theta), 5) + 105 * pow(cos(theta), 3) - 7 * cos(theta));
}


static int test_projection(void)
{
	struct sh_series *test = sh_series_new(20, 0);
	struct sh_series *exact = sh_series_new(20, 0);

	sh_series_from_func(test, Y_00_00, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 0, 0, 1);
	fprintf(stderr, "Projection of Y_{0,0} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_from_func(test, Y_07_01, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 7, 1, 1);
	fprintf(stderr, "Projection of Y_{7,1} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_from_func(test, Y_08_n7, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 8, -7, 1);
	fprintf(stderr, "Projection of Y_{8,-7} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_from_func(test, Y_10_03, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 10, 3, 1);
	fprintf(stderr, "Projection of Y_{10,3} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_free(test);
	sh_series_free(exact);

	return 0;
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


static void test1(int td, int fd, int speed, gsl_rng *rng, struct instrument_array *instruments, double delta_t, double *gmst, int n_gmst, int trials)
{
	struct correlator_network_baselines *baselines = correlator_network_baselines_new(instruments);

	const int sky_l_max = correlator_network_l_max(baselines, delta_t);
	const double dump_interval = correlator_dump_interval(sky_l_max, 20);
	const int time_series_length = dump_interval / delta_t;

	double injection[time_series_length];
	double *tseries[] = {
		malloc(time_series_length * sizeof(**tseries)),
		malloc(time_series_length * sizeof(**tseries))
	};
	complex double *fseries[] = {
		malloc(time_series_length * sizeof(**fseries)),
		malloc(time_series_length * sizeof(**fseries))
	};
	fftw_plan fftplans[] = {
		correlator_tseries_to_fseries_plan(tseries[0], fseries[0], time_series_length),
		correlator_tseries_to_fseries_plan(tseries[1], fseries[1], time_series_length)
	};

	struct correlator_network_plan_fd *fdplans = fd ? correlator_network_plan_fd_new(baselines, time_series_length, delta_t) : NULL;
	struct correlator_network_plan_td *tdplans = td ? correlator_network_plan_td_new(baselines, delta_t) : NULL;

	double *windows[] = {
		td ? correlator_square_window_new(time_series_length - 2 * tdplans->plans[0]->transient, 0, 1.0 / (time_series_length - 2 * tdplans->plans[0]->transient)) : NULL
	};

	struct sh_series *solution, *solution_lo;

	struct sh_series *tdsky = sh_series_new_zero(sky_l_max, 0);
	struct sh_series *fdsky = sh_series_new_zero(sky_l_max, 0);
	struct sh_series *tdaverage = sh_series_new(sky_l_max, 0);
	struct sh_series *fdaverage = sh_series_new(sky_l_max, 0);

	const double variance = 1.0;

	struct timeval tstart, tend;
	double runtime;

	int N_T, l_T, l_xi;

	int i, j, k;

	if(instrument_array_len(instruments) != 2) {
		fprintf(stderr, "FIXME!!\n");
		return;
	}

	fprintf(stderr, "=== Parameters ===\n");
	if(td) {
		fprintf(stderr, "=== Time-Domain Correlator ===\n");
		diagnostics_dump_correlator_plan_td_stats(stderr, tdplans->plans[0]);
	}
	if(fd) {
		fprintf(stderr, "=== Frequency-Domain Correlator ===\n");
		diagnostics_dump_correlator_plan_fd_stats(stderr, fdplans->plans[0]);
	}
	fprintf(stderr, "=== Other Numbers ===\n");
	fprintf(stderr, "l_{xi} = %d\n", sky_l_max);
	fprintf(stderr, "=== End Parameters ===\n");

	/* extract correlator parameters */

	N_T = tdplans->plans[0]->proj_a->n;
	l_T = tdplans->plans[0]->proj_a->l_max;
	l_xi = sky_l_max;

	/* loop over GMST */
	for(i = 0; i < n_gmst; i++) {
		fprintf(stderr, "=== Begin Test ===\nGMST = %g\n", gmst[i]);

		/* zero averages */
		sh_series_zero(tdaverage);
		sh_series_zero(fdaverage);

		/* record the start time */
		gettimeofday(&tstart, NULL);

		/* repeated trials */
		for(j = 0; j < trials; j++) {
			fprintf(stderr, "trial %d / %d\n", j + 1, trials);

			/* zero the instrument time series */
			for(k = 0; k < instrument_array_len(instruments); k++)
				memset(tseries[k], 0, time_series_length * sizeof(*tseries[k]));

			/* inject point source at vernal equinox */
			if(!speed) {
				gaussian_white_noise(injection, time_series_length, variance, rng);
				for(k = 0; k < instrument_array_len(instruments); k++)
					inject_into(instrument_array_get(instruments, k), tseries[k], injection, time_series_length, delta_t, 0.0, 0.0, gmst[i]);
			}

			/* compute integrated cross power */

			/* time domain */
			if(td) {
				correlator_network_integrate_power_td(tdsky, tseries, time_series_length, windows, tdplans);
				sh_series_rotate_z(tdsky, tdsky, gmst[i]);
				sh_series_scale(tdaverage, (double) j / (j + 1));
				sh_series_add(tdaverage, 1.0 / (j + 1), tdsky);
			}

			/* frequency domain */
			if(fd) {
				for(k = 0; k < instrument_array_len(instruments); k++)
					correlator_tseries_to_fseries(tseries[k], fseries[k], time_series_length, fftplans[k]);
				correlator_network_integrate_power_fd(fdsky, fseries, fdplans);
				sh_series_rotate_z(fdsky, fdsky, gmst[i]);
				sh_series_scale(fdaverage, (double) j / (j + 1));
				sh_series_add(fdaverage, 1.0 / (j + 1), fdsky);
			}
		}

		/* record the end time */
		gettimeofday(&tend, NULL);

		/* compute the execution speed */
		runtime = (tend.tv_sec - tstart.tv_sec) + 1e-6 * (tend.tv_usec - tstart.tv_usec);
		fprintf(stderr, "=== Begin Speed ===\n");
		fprintf(stderr, "%ld samples analyzed in %g seconds = %g samples/second\n", ((long int) time_series_length) * trials, runtime, ((long int) time_series_length) * trials / runtime);
		fprintf(stderr, "=== End Speed ===\n");

		/* compute "exact" solution for source at vernal equinox */
		solution = gaussian_noise_solution(baselines->baselines[0], 2 * sky_l_max, variance, delta_t, 0.0, 0.0, gmst[i]);

		/* compute exact solution for only those coefficients
		 * generated by the correlator */
		{
		int l, m;
		solution_lo = sh_series_copy(solution);
		for(l = sky_l_max + 1; l <= (int) solution_lo->l_max; l++)
			for(m = -l; m <= l; m++)
				sh_series_set(solution_lo, l, m, 0.0);
		}

		/* compute normalized RMS error */
		fprintf(stderr, "=== Begin Errors ===\n");
		{
		double td_fractional_rms_error;
		double fd_fractional_rms_error;
		double fractional_rms_error_lxi;
		double td_fractional_rms_error_N;
		double fd_fractional_rms_error_N;

		sh_series_clip(tdaverage, 1e-10);
		sh_series_clip(fdaverage, 1e-10);

		td_fractional_rms_error = diagnostics_rms_error(tdaverage, solution) / variance;
		fd_fractional_rms_error = diagnostics_rms_error(fdaverage, solution) / variance;
		fractional_rms_error_lxi = diagnostics_rms_error(solution_lo, solution) / variance;
		td_fractional_rms_error_N = diagnostics_rms_error(tdaverage, solution_lo) / variance;
		fd_fractional_rms_error_N = diagnostics_rms_error(fdaverage, solution_lo) / variance;

		/* error output */
		fprintf(stderr, "=== Time-Domain Correlator Errors ===\n");
		fprintf(stderr, "\\Delta \\bar{\\xi}_{\\mathrm{RMS}} = %g dB\n", dB(td_fractional_rms_error));
		fprintf(stderr, "\\Delta \\bar{\\xi}_{l_{\\xi} \\mathrm{RMS}} = %g dB\n", dB(fractional_rms_error_lxi));
		fprintf(stderr, "\\Delta \\bar{\\xi}_{\\delay \\mathrm{RMS}} = %g dB\n", dB(td_fractional_rms_error_N));
		fprintf(stderr, "td fractional error @ source = %g dB\n", dB((creal(sh_series_eval(solution, 0.0, 0.0)) - creal(sh_series_eval(tdaverage, 0.0, 0.0))) / variance));

		fprintf(stderr, "=== Frequency-Domain Correlator Errors ===\n");
		fprintf(stderr, "\\Delta \\bar{\\xi}_{\\mathrm{RMS}} = %g dB\n", dB(fd_fractional_rms_error));
		fprintf(stderr, "\\Delta \\bar{\\xi}_{l_{\\xi} \\mathrm{RMS}} = %g dB\n", dB(fractional_rms_error_lxi));
		fprintf(stderr, "\\Delta \\bar{\\xi}_{\\delay \\mathrm{RMS}} = %g dB\n", dB(fd_fractional_rms_error_N));
		fprintf(stderr, "fd fractional error @ source = %g dB\n", dB((creal(sh_series_eval(solution, 0.0, 0.0)) - creal(sh_series_eval(fdaverage, 0.0, 0.0))) / variance));
		fprintf(stderr, "fractional RMS (td - fd) = %g dB\n", dB(diagnostics_rms_error(fdaverage, tdaverage) / variance));
		}
		fprintf(stderr, "=== End Errors ===\n");

		/* solution output */
		{
		char filename[100];
		sprintf(filename, "tests_tdaverage_%d_%d_%d_%020.17f.dat", N_T, l_T, l_xi, gmst[i]);
		diagnostics_dump_sh_series(tdaverage, filename);
		sprintf(filename, "tests_fdaverage_%d_%d_%d_%020.17f.dat", N_T, l_T, l_xi, gmst[i]);
		diagnostics_dump_sh_series(fdaverage, filename);
		sprintf(filename, "tests_exact_%d_%d_%d_%020.17f.dat", N_T, l_T, l_xi, gmst[i]);
		diagnostics_dump_sh_series(solution, filename);
		}

		/* prepare for next test */
		sh_series_free(solution);
		sh_series_free(solution_lo);

		fprintf(stderr, "=== End Test ===\n");
	}

	/*
	 * Clean up
	 */

	for(k = 0; k < instrument_array_len(instruments); k++) {
		free(tseries[k]);
		free(fseries[k]);
		fftw_destroy_plan(fftplans[k]);
	}
	for(k = 0; k < baselines->n_baselines; k++)
		free(windows[k]);
	correlator_network_plan_fd_free(fdplans);
	correlator_network_plan_td_free(tdplans);
	correlator_network_baselines_free(baselines);
	sh_series_free(tdsky);
	sh_series_free(fdsky);
	sh_series_free(tdaverage);
	sh_series_free(fdaverage);
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
	instrument_array_append(instruments, instrument_new_from_r_theta_phi(+0.005, M_PI_2, 0, NULL, NULL));
	instrument_array_append(instruments, instrument_new_from_r_theta_phi(-0.005, M_PI_2, 0, NULL, NULL));

	const double delta_t = 1.0 / options->sample_frequency;

	const int n_gmst = 3;
	double gmst[n_gmst];
	int i, k;

	test_projection();

	gsl_rng_set(rng, time(NULL));

	for(i = 0; i < n_gmst; i++)
		gmst[i] = i * (M_PI / 4);

	test1(1, 1, 0, rng, instruments, delta_t, gmst, n_gmst, options->trials);


	/*
	 * Clean up
	 */

	instrument_array_free(instruments);
	command_line_options_free(options);
	gsl_rng_free(rng);

	return 0;
}
