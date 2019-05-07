/*
 * Copyright (C) 2007--2009,2012,2019  Kipp Cannon
 * Copyright (C) 2007  Patrick Sutton
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
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/inject.h>
#include <sphradiometer/correlator.h>

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


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


struct options {
	struct instrument_array *instruments;
	char *data1_cache_name;
	char *data1_channel_name;
	char *mdc1_cache_name;
	char *mdc1_channel_name;
	struct instrument *inst2;
	char *data2_cache_name;
	char *data2_channel_name;
	char *mdc2_cache_name;
	char *mdc2_channel_name;
	LIGOTimeGPS analysis_start;
	double analysis_duration;
	double integration_duration;
	char *output_name;
	double pass_band_low;	/* Hz, negative == disable */
	double pass_band_high;	/* Hz, negative == disable */
	double down_sample_to;	/* Hz, negative == disable */
};


static struct options *command_line_options_new(void)
{
	struct options *options = malloc(sizeof(*options));


	/*
	 * Defaults.
	 */


	*options = (struct options) {
		.instruments = instrument_array_new(2),
		.data1_cache_name = NULL,
		.data1_channel_name = NULL,
		.mdc1_cache_name = NULL,
		.mdc1_channel_name = NULL,
		.data2_cache_name = NULL,
		.data2_channel_name = NULL,
		.mdc2_cache_name = NULL,
		.mdc2_channel_name = NULL,
		.analysis_start = (LIGOTimeGPS) {
			.gpsSeconds = 0,
			.gpsNanoSeconds = 0,
		},
		.analysis_duration = 0.0,
		.integration_duration = 0.0,
		.output_name = NULL,
		.pass_band_low = -1.0,
		.pass_band_high = -1.0,
		.down_sample_to = -1.0,
	};

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
		instrument_array_free(options->instruments);
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
		{"data1-cache-name",		required_argument,	NULL,	'A'},
		{"data1-channel-name",		required_argument,	NULL,	'B'},
		{"mdc1-cache-name",		required_argument,	NULL,	'C'},
		{"mdc1-channel-name",		required_argument,	NULL,	'D'},
		{"data2-cache-name",		required_argument,	NULL,	'a'},
		{"data2-channel-name",		required_argument,	NULL,	'b'},
		{"mdc2-cache-name",		required_argument,	NULL,	'c'},
		{"mdc2-channel-name",		required_argument,	NULL,	'd'},
		{"analysis-start-gps",		required_argument,	NULL,	'E'},
		{"analysis-duration",		required_argument,	NULL,	'F'},
		{"integration-duration",	required_argument,	NULL,	'G'},
		{"output",			required_argument,	NULL,	'H'},
		{"pass-band-low",		required_argument,	NULL,	'I'},
		{"pass-band-high",		required_argument,	NULL,	'J'},
		{"down-sample-to",		required_argument,	NULL,	'K'},
		{"help",			no_argument,		NULL,	'h'},
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
	case 'B': {
		char inst_name[3] = {optarg[0], optarg[1], '\0'};
		instrument_array_set(options->instruments, 0, instrument_new_from_name(inst_name));
		options->data1_channel_name = optarg;
		break;
	}

	/* mdc1-cache-name */
	case 'C':
		options->mdc1_cache_name = optarg;
		break;

	/* mdc1-channel-name */
	case 'D':
		options->mdc1_channel_name = optarg;
		break;

	/* data2-cache-name */
	case 'a':
		options->data2_cache_name = optarg;
		break;

	/* data2-channel-name */
	case 'b': {
		char inst_name[3] = {optarg[0], optarg[1], '\0'};
		instrument_array_set(options->instruments, 1, instrument_new_from_name(inst_name));
		options->data2_channel_name = optarg;
		break;
	}

	/* mdc2-cache-name */
	case 'c':
		options->mdc2_cache_name = optarg;
		break;

	/* mdc2-channel-name */
	case 'd':
		options->mdc2_channel_name = optarg;
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

	/* output */
	case 'H':
		options->output_name = optarg;
		break;

	/* pass band low frequency limit */
	case 'I':
		options->pass_band_low = atof(optarg);
		break;

	/* pass band high frequency limit */
	case 'J':
		options->pass_band_high = atof(optarg);
		break;

	/* down sampling target frequency */
	case 'K':
		options->down_sample_to = atof(optarg);
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


static REAL8TimeSeries *get_real8series_from_cache(
	const char *cache_name,
	const char *channel_name,
	LIGOTimeGPS start,
	double duration,
	size_t lengthlimit
)
{
	LALCache *cache;
	LALFrStream *stream;
	REAL8TimeSeries *data;
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
	data = XLALFrStreamReadREAL8TimeSeries(stream, channel_name, &start, duration, lengthlimit);

	/* check for gaps and close */
	gap = stream->state & LAL_FR_STREAM_GAP;
	XLALFrStreamClose(stream);

	/* error checking */
	if(!data)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if(gap) {
		XLALDestroyREAL8TimeSeries(data);
		XLALPrintError("error: gap detected in input data");
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/* done */
	return data;
}


static REAL8TimeSeries *get_instrument_time_series(
	const char *data_cache_name,
	const char *data_channel_name,
	const char *mdc_cache_name,
	const char *mdc_channel_name,
	LIGOTimeGPS start,
	double duration,
	size_t lengthlimit
)
{
	REAL8TimeSeries *series;

	/* get data */
	series = get_real8series_from_cache(data_cache_name, data_channel_name, start, duration, lengthlimit);
	if(!series)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* add mdc injection */
	if(mdc_cache_name && mdc_channel_name) {
		REAL8TimeSeries *mdc;
		int i;

		mdc = get_real8series_from_cache(mdc_cache_name, mdc_channel_name, start, duration, lengthlimit);
		if(!mdc) {
			XLALDestroyREAL8TimeSeries(series);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		for(i = 0; i < (int) series->data->length; i++)
			series->data->data[i] += mdc->data->data[i];

		XLALDestroyREAL8TimeSeries(mdc);
	}

	return series;
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
 *                          Network Response Matrix
 *
 * ============================================================================
 */


static void FplusDP(double *fplus1, double *fplus2, const LALDetector *det1, const LALDetector *det2, double theta, double phi)
{
	double fcross1, fcross2;
	double psi;

	XLALComputeDetAMResponse(fplus1, &fcross1, det1->response, phi, M_PI_2 - theta, 0.0, 0.0);
	XLALComputeDetAMResponse(fplus2, &fcross2, det2->response, phi, M_PI_2 - theta, 0.0, 0.0);

	psi = atan(2.0 * (*fplus1 * fcross1 + *fplus2 * fcross2) / (*fplus1 * *fplus1 + *fplus2 * *fplus2 - fcross1 * fcross1 - fcross2 * fcross2)) / 4;

	*fplus1 = cos(2 * psi) * *fplus1 + sin(2 * psi) * fcross1;
	*fplus2 = cos(2 * psi) * *fplus2 + sin(2 * psi) * fcross2;
}


struct response_element_data {
	const LALDetector *det1, *det2;
};


static double response12_element(double theta, double phi, void *data)
{
	double fplus1, fplus2;

	FplusDP(&fplus1, &fplus2, ((struct response_element_data *) data)->det1, ((struct response_element_data *) data)->det2, theta, phi);

	return fplus1 * fplus2 / (fplus1 * fplus1 + fplus2 * fplus2);
}


/* FIXME: We should be picking up the instruments from a baseline structure so
 * that we know we've got the right ones */
static struct sh_series *response12(const LALDetector *det1, const LALDetector *det2, unsigned int l_max)
{
	struct sh_series *series = sh_series_new(l_max, 0);
	struct response_element_data data = {
		.det1 = det1,
		.det2 = det2,
	};

	return sh_series_from_realfunc(series, response12_element, &data);
}


static struct correlator_plan_fd *correlator_plan_mult_by_response(struct correlator_plan_fd *plan, const struct sh_series *response)
{
	struct sh_series *result;
	struct sh_series_product_plan *product_plan;
	int i;

	result = sh_series_new(plan->delay_product->l_max, 0);
	if(!result)
		return NULL;
	product_plan = sh_series_product_plan_new(result, &plan->delay_product->series[0], response);
	if(!product_plan) {
		sh_series_free(result);
		return NULL;
	}

	for(i = 0; i < plan->delay_product->n; i++) {
		sh_series_product(result, &plan->delay_product->series[i], response, product_plan);
		sh_series_assign(&plan->delay_product->series[i], result);
	}

	sh_series_product_plan_free(product_plan);
	sh_series_free(result);

	return plan;
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
	REAL8TimeSeries *series[2];
	struct correlator_network_baselines *baselines;
	int time_series_length;
	double *tseries;
	complex double *fseries[2];
	fftw_plan fftplans[2];
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


	series[0] = get_instrument_time_series(options->data1_cache_name, options->data1_channel_name, options->mdc1_cache_name, options->mdc1_channel_name, options->analysis_start, options->analysis_duration, 0);
	series[1] = get_instrument_time_series(options->data2_cache_name, options->data2_channel_name, options->mdc2_cache_name, options->mdc2_channel_name, options->analysis_start, options->analysis_duration, 0);

	if(!series[0] || !series[1]) {
		XLALPrintError("failure loading data\n");
		exit(1);
	}
	if(series[0]->deltaT != series[1]->deltaT) {
		fprintf(stderr, "sample frequency mismatch (%g != %g)\n", series[0]->deltaT, series[1]->deltaT);
		exit(1);
	}

	series[0] = condition_time_series(series[0], options->pass_band_low, options->pass_band_high, options->down_sample_to);
	series[1] = condition_time_series(series[1], options->pass_band_low, options->pass_band_high, options->down_sample_to);

	time_series_length = options->integration_duration / series[0]->deltaT;
	if(time_series_length * series[0]->deltaT != options->integration_duration) {
		fprintf(stderr, "integration duration not an integer number of samples\n");
		exit(1);
	}

	tseries = malloc(time_series_length * sizeof(*tseries));
	fseries[0] = malloc(time_series_length * sizeof(**fseries));
	fseries[1] = malloc(time_series_length * sizeof(**fseries));

	fftplans[0] = correlator_tseries_to_fseries_plan(tseries, fseries[0], time_series_length);
	fftplans[1] = correlator_tseries_to_fseries_plan(tseries, fseries[1], time_series_length);


	/*
	 * Construct the correlator.
	 */


	baselines = correlator_network_baselines_new(options->instruments);
	fdplans = correlator_network_plan_fd_new(baselines, time_series_length, series[0]->deltaT);
	sky_l_max = correlator_network_l_max(baselines, series[0]->deltaT);
	sky = sh_series_new_zero(sky_l_max, 0);


	/*
	 * Shoe-horn the response matrix into the correlator.
	 * FIXME: this is ugly
	 */


	{
	/* turn off azimuthal symmetry flag */
	struct sh_series_array *new = sh_series_array_new(fdplans->plans[0]->delay_product->n, fdplans->plans[0]->delay_product->l_max, 0);
	sh_series_array_assign(new, fdplans->plans[0]->delay_product);
	sh_series_array_free(fdplans->plans[0]->delay_product);
	fdplans->plans[0]->delay_product = new;
	}
	{
	/* construct baseline response function */
	struct sh_series *resp12 = response12(instrument_array_get(options->instruments, 0)->data, instrument_array_get(options->instruments, 1)->data, 10);
	/* multiply each element of the correlator's transform matrix by
	 * the response */
	correlator_plan_mult_by_response(fdplans->plans[0], resp12);
	/* clean up */
	sh_series_free(resp12);
	}


	/*
	 * Loop over integration chunk.
	 */


	fprintf(stderr, "starting integration\n");
	gettimeofday(&t_start, NULL);
	for(start_sample = 0; start_sample + time_series_length <= (int) series[0]->data->length; start_sample += time_series_length / 2) {
		/*
		 * Extract the data to integrate.
		 */


		fprintf(stderr, "\t[%.3f s, %.3f s)\n", start_sample * series[0]->deltaT, (start_sample + time_series_length) * series[0]->deltaT);
		for(k = 0; k < instrument_array_len(options->instruments); k++) {
			memcpy(tseries, &series[k]->data->data[start_sample], time_series_length * sizeof(*tseries));
			correlator_tseries_to_fseries(tseries, fseries[k], time_series_length, fftplans[k]);
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


		/* FIXME: should probably do something here */
	}
	gettimeofday(&t_end, NULL);
	fprintf(stderr, "finished integration\n");

	fprintf(stderr, "analyzed %g s of data in %g s chunks, overlapping 50%%, in %g s\n", series[0]->data->length * series[0]->deltaT, time_series_length * series[0]->deltaT, (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec) * 1e-6);
	fprintf(stderr, "data sample rate was %g Hz\n", 1.0 / series[0]->deltaT);
	fprintf(stderr, "sky was computed to harmonic order l = %d\n", sky->l_max);


	/*
	 * Clean up
	 */


	free(tseries);
	for(k = 0; k < instrument_array_len(options->instruments); k++) {
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
