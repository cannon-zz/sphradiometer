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
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <radiometer/instrument.h>
#include <radiometer/sh_series.h>
#include <radiometer/inject.h>
#include <radiometer/correlator.h>
#include <diagnostics.h>
#include <instruments.h>


/*
 * ============================================================================
 *
 *                                Diagnostics
 *
 * ============================================================================
 */


static void dump_network_stats(FILE *f, const struct correlator_network_baselines *baselines, const struct correlator_network_plan_td *tdplans, const struct correlator_network_plan_fd *fdplans)
{
	int i;

	for(i = 0; i < baselines->n_baselines; i++) {
		fprintf(f, "=== Baseline %d ===\n", i);
		fprintf(f, "Time Domain Correlator:n");
		diagnostics_dump_correlator_plan_td_stats(f, tdplans->plans[i]);
		fprintf(f, "Frequency Domain Correlator:n");
		diagnostics_dump_correlator_plan_fd_stats(f, fdplans->plans[i]);
	}
	fprintf(f, "=== End Baselines ===\n");
}


static void dump_time_series(const double *series, int n, char *name)
{
	FILE *f = fopen(name, "w");
	int i;

	for(i = 0; i < n; i++)
		fprintf(f, "%d\t%g\n", i, *series++);

	fclose(f);
}


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


struct options {
	struct {
		double ra;
		double dec;
	} *injection_ra_dec;
	int n_injections;
	double baseline_theta;
	double baseline_phi;
};


static struct options *command_line_options_new(void)
{
	struct options *options = malloc(sizeof(*options));

	/* defaults */
	*options = (struct options) {
		.injection_ra_dec = NULL,
		.n_injections = 0,
		.baseline_theta = M_PI_2,
		.baseline_phi = 0,
	};

	return options;
}


static void command_line_options_free(struct options *options)
{
	if(options)
		free(options->injection_ra_dec);
	free(options);
}


static int command_line_add_injection(struct options *options, double ra, double dec)
{
	void *new = realloc(options->injection_ra_dec, ++options->n_injections * sizeof(*options->injection_ra_dec));

	if(!new)
		return -1;

	options->injection_ra_dec = new;
	options->injection_ra_dec[options->n_injections - 1].ra = ra;
	options->injection_ra_dec[options->n_injections - 1].dec = dec;

	return 0;
}


struct options *parse_command_line(int argc, char *argv[])
{
	int c;
	int option_index;
	struct option long_options[] = {
		{"injection-ra-dec",	required_argument,	NULL,	'A'},
		{"baseline-theta",	required_argument,	NULL,	'B'},
		{"baseline-phi",	required_argument,	NULL,	'C'},
		{NULL,	0,	NULL,	0}
	};
	struct options *options = command_line_options_new();

	if(!options)
		return NULL;

	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
	/* injection-ra-dec */
	case 'A': {
		double ra;
		double dec;
		sscanf(optarg, "%lf,%lf", &ra, &dec);
		command_line_add_injection(options, ra, dec);
		break;
	}

	/* baseline-theta */
	case 'B':
		options->baseline_theta = atof(optarg);
		break;

	/* baseline-phi */
	case 'C':
		options->baseline_phi = atof(optarg);
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
 *                                Entry Point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	struct options *options = parse_command_line(argc, argv);
	const double delta_t = 1.0 / 512;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);

	struct instrument *instruments[] = {
		instrument_from_name("H1"),
		instrument_from_name("L1")/*,
		instrument_from_name("G1")*/
	};
	int n_instruments = sizeof(instruments) / sizeof(*instruments);

	struct correlator_network_baselines *baselines = correlator_network_baselines_new(instruments, n_instruments);

	const int sky_l_max = correlator_network_l_max(baselines, delta_t);
	const double dump_interval = correlator_dump_interval(sky_l_max, 20);
	const int time_series_length = dump_interval / delta_t;

	double injection[time_series_length];
	double *tseries[] = {
		malloc(time_series_length * sizeof(**tseries)),
		malloc(time_series_length * sizeof(**tseries)),
		malloc(time_series_length * sizeof(**tseries))
	};
	complex double *fseries[] = {
		malloc(time_series_length * sizeof(**fseries)),
		malloc(time_series_length * sizeof(**fseries)),
		malloc(time_series_length * sizeof(**fseries))
	};
	fftw_plan fftplans[] = {
		correlator_tseries_to_fseries_plan(tseries[0], fseries[0], time_series_length),
		correlator_tseries_to_fseries_plan(tseries[1], fseries[1], time_series_length),
		correlator_tseries_to_fseries_plan(tseries[2], fseries[2], time_series_length)
	};

	struct correlator_network_plan_fd *fdplans = correlator_network_plan_fd_new(baselines, time_series_length, delta_t);
	struct correlator_network_plan_td *tdplans = correlator_network_plan_td_new(baselines, delta_t);

	double *windows[] = {
		correlator_square_window_new(time_series_length - 2 * tdplans->plans[0]->transient, 0, 1.0 / (time_series_length - 2 * tdplans->plans[0]->transient))/*,
		correlator_square_window_new(time_series_length - 2 * tdplans->plans[1]->transient, 0, 1.0 / (time_series_length - 2 * tdplans->plans[1]->transient)),
		correlator_square_window_new(time_series_length - 2 * tdplans->plans[2]->transient, 0, 1.0 / (time_series_length - 2 * tdplans->plans[2]->transient))*/
	};

	struct sh_series *tdsky = sh_series_new_zero(sky_l_max, 0);
	struct sh_series *fdsky = sh_series_new_zero(sky_l_max, 0);
	struct sh_series *tdaverage = sh_series_new_zero(sky_l_max, 0);
	struct sh_series *fdaverage = sh_series_new_zero(sky_l_max, 0);
	const int n_dumps = 86164 / dump_interval;
	const double delta_gmst = 2 * M_PI / n_dumps;
	double gmst;
	char filename[100];
	int i, j, k;

	dump_network_stats(stderr, baselines, tdplans, fdplans);
	fprintf(stderr, "sky: l max = %d\n", sky_l_max);
	fprintf(stderr, "data: time series length = %g s (%d samples)\n", dump_interval, time_series_length);

	gsl_rng_set(rng, time(NULL));

	/* loop over GMST */
	for(i = 0, gmst = 0.0; i < n_dumps; gmst = ++i * delta_gmst) {
		fprintf(stderr, "GMST = %g rad\n", gmst);

		/* zero the instrument time series */
		for(k = 0; k < n_instruments; k++)
			memset(tseries[k], 0, time_series_length * sizeof(*tseries[k]));

		/* generate instrumental noise for 1e-3 SNR injections *
		gaussian_white_noise(tseries_a, time_series_length, 1000.0, rng);
		gaussian_white_noise(tseries_b, time_series_length, 1000.0, rng);*/

		/* inject point sources */
		for(j = 0; j < options->n_injections; j++) {
			gaussian_white_noise(injection, time_series_length, 1.0, rng);
			for(k = 0; k < n_instruments; k++)
				inject_into(instruments[k], tseries[k], injection, time_series_length, delta_t, options->injection_ra_dec[j].ra, options->injection_ra_dec[j].dec, gmst);
		}

		/* compute integrated cross power */
		/* time domain */
		/*correlator_network_integrate_power_td(tdsky, tseries, time_series_length, windows, tdplans);
		sh_series_rotate_z(tdsky, tdsky, gmst);*/

		/* frequency domain */
		for(k = 0; k < n_instruments; k++)
			correlator_tseries_to_fseries(tseries[k], fseries[k], time_series_length, fftplans[k]);
		correlator_network_integrate_power_fd(fdsky, fseries, fdplans);
		sh_series_rotate_z(fdsky, fdsky, gmst);

		/* output */
		sprintf(filename, "power_td_%020.17f.dat", gmst);
		diagnostics_dump_sh_series(tdsky, filename);
		sprintf(filename, "power_fd_%020.17f.dat", gmst);
		diagnostics_dump_sh_series(fdsky, filename);

		sh_series_scale(tdaverage, (double) i / (i + 1));
		sh_series_add(tdaverage, 1.0 / (i + 1), tdsky);
		sh_series_scale(fdaverage, (double) i / (i + 1));
		sh_series_add(fdaverage, 1.0 / (i + 1), fdsky);

		sprintf(filename, "averg_td_%020.17f.dat", gmst);
		diagnostics_dump_sh_series(tdaverage, filename);
		sprintf(filename, "averg_fd_%020.17f.dat", gmst);
		diagnostics_dump_sh_series(fdaverage, filename);
	}

	/*
	 * Clean up
	 */

	for(k = 0; k < n_instruments; k++) {
		instrument_free(instruments[k]);
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
	command_line_options_free(options);
	gsl_rng_free(rng);

	return 0;
}
