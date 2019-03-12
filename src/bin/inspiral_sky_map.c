/*
 * Copyright (C) 2009  Kipp Cannon, Chad Hanna
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
#include <sys/time.h>
#include <time.h>
#include <getopt.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <radiometer/instrument.h>
#include <radiometer/sh_series.h>
#include <radiometer/inject.h>
#include <radiometer/correlator.h>
#include <diagnostics.h>
#include <instruments.h>

#include <lal/BandPassTimeSeries.h>
#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimulation.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/LALError.h>
#include <lal/XLALError.h>
#include <lal/Units.h>
#include <lal/Window.h>


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


struct options {
	struct instrument *instruments[3];
	const LALDetector *detectors[3];	/* FIXME: merge with instruments[] */
	int n_instruments;
	char *data1_cache_name;
	char *data1_channel_name;
	char *data2_cache_name;
	char *data2_channel_name;
	char *data3_cache_name;
	char *data3_channel_name;
	LIGOTimeGPS analysis_start;
	double analysis_duration;
	double integration_duration;
	double integration_shift;
	char *output_name;
};


static struct options *command_line_options_new(void)
{
	struct options *options = malloc(sizeof(*options));


	/*
	 * Defaults.
	 */


	*options = (struct options) {
		.instruments = {
			NULL,
			NULL,
			NULL,
		},
		.detectors = {
			NULL,
			NULL,
			NULL,
		},
		.n_instruments = 3,
		.data1_cache_name = NULL,
		.data1_channel_name = NULL,
		.data2_cache_name = NULL,
		.data2_channel_name = NULL,
		.analysis_start = (LIGOTimeGPS) {
			.gpsSeconds = 0,
			.gpsNanoSeconds = 0,
		},
		.analysis_duration = 0.0,
		.integration_duration = 0.0,
		.integration_shift = 0.125,
		.output_name = NULL,
	};

	return options;
}


static struct options *command_line_set_instrument(struct options *options, int n, const char *name)
{
	char instrument_name[3] = {name[0], name[1], '\0'};

	/* correct VIRGO's name appearing in CBC channel names */
	if(!strcmp(instrument_name, "V1"))
		instrument_name[1] = '2';

	options->detectors[n] = XLALDetectorPrefixToLALDetector(instrument_name);
	options->instruments[n] = instrument_from_LALDetector(options->detectors[n]);

	return options;
}


static int command_line_options_validate(struct options *options)
{
	/*
	 * All OK.
	 */

	return 0;
}


static void command_line_options_free(struct options *options)
{
	if(options) {
		int i;
		for(i = 0; i < options->n_instruments; i++)
			instrument_free(options->instruments[i]);
	}
	free(options);
}


static void usage(int argc, char *argv[])
{
	fprintf(stderr, "%s: no help message available at this time ...\n", argv[0]);
}


struct options *command_line_parse(int argc, char *argv[])
{
	int c;
	int option_index;
	struct option long_options[] = {
		{"data1-cache-name",	required_argument,	NULL,	'A'},
		{"data1-channel-name",	required_argument,	NULL,	'B'},
		{"data2-cache-name",	required_argument,	NULL,	'a'},
		{"data2-channel-name",	required_argument,	NULL,	'b'},
		{"data3-cache-name",	required_argument,	NULL,	'C'},
		{"data3-channel-name",	required_argument,	NULL,	'c'},
		{"analysis-start-gps",	required_argument,	NULL,	'E'},
		{"analysis-duration",	required_argument,	NULL,	'F'},
		{"integration-duration",	required_argument,	NULL,	'G'},
		{"integration-shift",	required_argument,	NULL,	'g'},
		{"output",		required_argument,	NULL,	'H'},
		{"help",		no_argument,		NULL,	'h'},
		{NULL,	0,	NULL,	0}
	};
	struct options *options = command_line_options_new();

	if(!options)
		return NULL;

	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
	/* data1-cache-name */
	case 'A':
		options->data1_cache_name = optarg;
		break;

	/* data1-channel-name */
	case 'B':
		command_line_set_instrument(options, 0, optarg);
		options->data1_channel_name = optarg;
		break;

	/* data2-cache-name */
	case 'a':
		options->data2_cache_name = optarg;
		break;

	/* data2-channel-name */
	case 'b':
		command_line_set_instrument(options, 1, optarg);
		options->data2_channel_name = optarg;
		break;

	/* data3-cache-name */
	case 'C':
		options->data3_cache_name = optarg;
		break;

	/* data3-channel-name */
	case 'c':
		command_line_set_instrument(options, 2, optarg);
		options->data3_channel_name = optarg;
		break;

	/* analysis-start-gps */
	case 'E':
		XLALStrToGPS(&options->analysis_start, optarg, NULL);
		break;

	/* analysis-duration */
	case 'F':
		options->analysis_duration = atof(optarg);
		break;

	/* integration-duration */
	case 'G':
		options->integration_duration = atof(optarg);
		break;

	/* integration-shift */
	case 'g':
		options->integration_shift = atof(optarg);
		break;

	/* output */
	case 'H':
		options->output_name = optarg;
		break;

	/* help */
	case 'h':
		usage(argc, argv);
		exit(0);

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

	if(command_line_options_validate(options)) {
		command_line_options_free(options);
		return NULL;
	}

	return options;
}


/*
 * ============================================================================
 *
 *                                 Data Input
 *
 * ============================================================================
 */


static REAL8TimeSeries *get_real8series_from_txt(
	const char *filename,
	const char *channel_name,
	LIGOTimeGPS start,
	double duration
)
{
	FILE *f;
	REAL8TimeSeries *series;
	int i;
	int column;

	if(!strncmp(channel_name, "H1", 2)) {
		column = 1;
	} else if(!strncmp(channel_name, "L1", 2)) {
		column = 3;
	} else if(!strncmp(channel_name, "V1", 2)) {
		column = 5;
	} else
		return NULL;

	f = fopen(filename, "r");

	series = XLALCreateREAL8TimeSeries(channel_name, &start, 0.0, 1.0 / 2048, &lalDimensionlessUnit, 1.0 * 2048);

	for(i = 0; i < series->data->length; i++) {
		double values[7];
		fscanf(f, "%lg %lg %lg %lg %lg %lg %lg", &values[0], &values[1], &values[2], &values[3], &values[4], &values[5], &values[6]);

		series->data->data[i] = cabs(values[column] + I * values[column + 1]);
	}

	fclose(f);

	return series;
}


static REAL8TimeSeries *get_real8series_from_cache(
	const char *cache_name,
	const char *channel_name,
	LIGOTimeGPS start,
	double duration
)
{
	LALCache *cache;
	LALFrStream *stream;
	REAL8TimeSeries *data;
	COMPLEX8TimeSeries *indata;
	unsigned i;
	int gap;

	/* construct stream */
	cache = XLALCacheImport(cache_name);
	if(!cache)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	stream = XLALFrStreamCacheOpen(cache);
	XLALDestroyCache(cache);
	if(!stream)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* turn on checking for missing data */
	stream->mode = LAL_FR_STREAM_VERBOSE_MODE;

	/* get data */
#if 0
	indata = XLALFrReadCOMPLEX8TimeSeries(stream, channel_name, &start, duration, 0);
#else
	indata = XLALCreateCOMPLEX8TimeSeries(channel_name, &start, 0.0, 1.0 / 4096, &lalDimensionlessUnit, 1.0 * 4096);
	XLALFrStreamGetCOMPLEX8TimeSeries(indata, stream);
	indata->epoch = start;
#endif

	/* check for gaps and close */
	gap = stream->state & LAL_FR_STREAM_GAP;
	XLALFrStreamClose(stream);

	/* error checking */
	if(!indata)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if(gap) {
		XLALDestroyCOMPLEX8TimeSeries(indata);
		XLALPrintError("error: gap detected in input data");
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/* copy magnitude time series */
	data = XLALCreateREAL8TimeSeries(indata->name, &indata->epoch, indata->f0, indata->deltaT, &indata->sampleUnits, indata->data->length);
	if(!data) {
		XLALDestroyCOMPLEX8TimeSeries(indata);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	for(i = 0; i < indata->data->length; i++)
		data->data->data[i] = cabsf(indata->data->data[i]);
	XLALDestroyCOMPLEX8TimeSeries(indata);

	/* done */
	return data;
}


/*
 * ============================================================================
 *
 *                             Data Conditioning
 *
 * ============================================================================
 */


static REAL8TimeSeries *condition_time_series(
	REAL8TimeSeries *series,
	double flow,
	double fhigh,
	double sample_rate
)
{
	double new_deltaT = sample_rate > 0 ? 1.0 / sample_rate : series->deltaT;

	if(flow > 0)
		if(XLALHighPassREAL8TimeSeries(series, flow, .9, 6))
			XLAL_ERROR_NULL(XLAL_EFUNC);

	if(fhigh > 0)
		if(XLALLowPassREAL8TimeSeries(series, fhigh, .9, 6))
			XLAL_ERROR_NULL(XLAL_EFUNC);

	if(fabs(new_deltaT - series->deltaT) / series->deltaT >= 1e-2)
		if(XLALResampleREAL8TimeSeries(series, new_deltaT))
			XLAL_ERROR_NULL(XLAL_EFUNC);

	return series;
}


/*
 * ============================================================================
 *
 *                              Time Conversion
 *
 * ============================================================================
 */


static double gmst_from_epoch_and_offset(LIGOTimeGPS epoch, double offset)
{
	return XLALGreenwichSiderealTime(XLALGPSAdd(&epoch, offset), 0.0);
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
	gsl_rng *rng;
	struct options *options;
	REAL8TimeSeries *series[3];
	struct correlator_network_baselines *baselines;
	double *tseries;
	REAL8Window *window;
	complex double *fseries[3];
	fftw_plan fftplans[3];
	struct correlator_network_plan_fd *fdplans;
	int sky_l_max;
	struct sh_series *sky;
	double gmst;
	int start_sample;
	int k;
	struct timeval t_start, t_end;


	/*
	 * Parse command line.
	 */


	options = command_line_parse(argc, argv);


	/*
	 * Initialize the random number generator.
	 */


	rng = gsl_rng_alloc(gsl_rng_ranlxd1);
	gsl_rng_set(rng, time(NULL));


	/*
	 * Load time series data.
	 */


/*
	series[0] = get_real8series_from_cache(options->data1_cache_name, options->data1_channel_name, options->analysis_start, options->analysis_duration);
	series[1] = get_real8series_from_cache(options->data2_cache_name, options->data2_channel_name, options->analysis_start, options->analysis_duration);
	series[2] = get_real8series_from_cache(options->data3_cache_name, options->data3_channel_name, options->analysis_start, options->analysis_duration);
*/
	series[0] = get_real8series_from_txt("TestInjection1.txt", options->data1_channel_name, options->analysis_start, options->analysis_duration);
	series[1] = get_real8series_from_txt("TestInjection1.txt", options->data2_channel_name, options->analysis_start, options->analysis_duration);
	series[2] = get_real8series_from_txt("TestInjection1.txt", options->data3_channel_name, options->analysis_start, options->analysis_duration);

	if(!series[0] || !series[1] || !series[2]) {
		XLALPrintError("failure loading data\n");
		exit(1);
	}
	if(series[0]->deltaT != series[1]->deltaT || series[0]->deltaT != series[2]->deltaT) {
		fprintf(stderr, "sample frequency mismatch\n");
		exit(1);
	}

/*
	series[0] = condition_time_series(series[0], -1, -1, 1024);
	series[1] = condition_time_series(series[1], -1, -1, 1024);
	series[2] = condition_time_series(series[2], -1, -1, 1024);
*/

	window = XLALCreateRectangularREAL8Window(options->integration_duration / series[0]->deltaT);
	if(window->data->length * series[0]->deltaT != options->integration_duration) {
		fprintf(stderr, "integration duration not an integer number of samples\n");
		exit(1);
	}

	tseries = malloc(window->data->length * sizeof(*tseries));
	fseries[0] = malloc(window->data->length * sizeof(**fseries));
	fseries[1] = malloc(window->data->length * sizeof(**fseries));
	fseries[2] = malloc(window->data->length * sizeof(**fseries));

	fftplans[0] = correlator_tseries_to_fseries_plan(tseries, fseries[0], window->data->length);
	fftplans[1] = correlator_tseries_to_fseries_plan(tseries, fseries[1], window->data->length);
	fftplans[2] = correlator_tseries_to_fseries_plan(tseries, fseries[2], window->data->length);


	/*
	 * Construct the correlator.
	 */


	baselines = correlator_network_baselines_new(options->instruments, options->n_instruments);
	fdplans = correlator_network_plan_fd_new(baselines, window->data->length, series[0]->deltaT);
	sky_l_max = correlator_network_l_max(baselines, series[0]->deltaT);
	sky = sh_series_new_zero(sky_l_max, 0);


	/*
	 * Loop over integration chunk.
	 */


	fprintf(stderr, "starting integration\n");
	gettimeofday(&t_start, NULL);
	for(start_sample = 0; start_sample + window->data->length <= series[0]->data->length; start_sample += window->data->length * options->integration_shift) {
		char filename[128];

		/*
		 * Extract the data to integrate.
		 */


		fprintf(stderr, "\t[%.3f s, %.3f s)\n", start_sample * series[0]->deltaT, (start_sample + window->data->length) * series[0]->deltaT);
		for(k = 0; k < options->n_instruments; k++) {
			unsigned j;
			memcpy(tseries, &series[k]->data->data[start_sample], window->data->length * sizeof(*tseries));
			for(j = 0; j < window->data->length; j++)
				tseries[j] *= window->data->data[j] * sqrt(window->data->length / window->sumofsquares);
			correlator_tseries_to_fseries(tseries, fseries[k], window->data->length, fftplans[k]);
		}


		/*
		 * Compute angular distribution of integrated cross power.
		 */


		correlator_network_integrate_power_fd(sky, fseries, fdplans);

		/*
		 * Rotate sky.
		 */


		gmst = gmst_from_epoch_and_offset(options->analysis_start, start_sample * series[0]->deltaT);
		sh_series_rotate_z(sky, sky, gmst);


		/*
		 * Output.
		 */

		sprintf(filename, "%s_%020.17f.dat", options->output_name, gmst);
		diagnostics_dump_sh_series(sky, filename);
	}
	gettimeofday(&t_end, NULL);
	fprintf(stderr, "finished integration\n");

	fprintf(stderr, "analyzed %g s of data in %g s chunks, overlapping %g%%, in %g s\n", series[0]->data->length * series[0]->deltaT, window->data->length * series[0]->deltaT, 100 * (1 - options->integration_shift), (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec) * 1e-6);
	fprintf(stderr, "data sample rate was %g Hz\n", 1.0 / series[0]->deltaT);
	fprintf(stderr, "sky was computed to harmonic order l = %d\n", sky->l_max);


	/*
	 * Clean up
	 */


	free(tseries);
	for(k = 0; k < options->n_instruments; k++) {
		free(fseries[k]);
		fftw_destroy_plan(fftplans[k]);
	}
	correlator_network_plan_fd_free(fdplans);
	correlator_network_baselines_free(baselines);
	sh_series_free(sky);
	command_line_options_free(options);
	gsl_rng_free(rng);

	return 0;
}
