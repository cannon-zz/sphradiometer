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
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/correlator.h>
#include <sphradiometer/diagnostics.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALCache.h>
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimulation.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/Units.h>
#include <lal/Window.h>


#define Projection_lmax 8


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


struct options {
	struct instrument_array *instruments;
	char *snr_cache;
	char *noise_cache;
	char **channels;
	char *output_name;
};


static struct options *command_line_options_new(void)
{
	struct options *options = malloc(sizeof(*options));


	/*
	 * Defaults.
	 */


	*options = (struct options) {
		.instruments = instrument_array_new(0),
		.snr_cache = NULL,
		.channels = NULL,
		.noise_cache = NULL,
		.output_name = NULL,
	};

	return options;
}


static struct options *command_line_set_instrument(struct options *options, const char *name)
{
	char instrument_name[3] = {name[0], name[1], '\0'};
	struct instrument *instrument = instrument_new_from_name(instrument_name);

	if(!instrument) {
		XLALPrintError("failed to initiaize for instrument \"%s\"", instrument_name);
		return NULL;
	}

	if(!instrument_array_append(options->instruments, instrument)) {
		XLALPrintError("failed to append instrument \"%s\" to list", instrument_name);
		return NULL;
	}

	options->channels = realloc(options->channels,  instrument_array_len(options->instruments) * sizeof(*options->channels));
	if(!options->channels) {
		/* FIXME:  leaks memory */
		XLALPrintError("failed to resize channel list");
		return NULL;
	}
	options->channels[instrument_array_len(options->instruments) - 1] = optarg;

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
		/* don't free optarg memory, just the array of pointers to them */
		free(options->channels);
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
		{"snr-cache",	required_argument,	NULL,	'A'},
		{"snr-channel",	required_argument,	NULL,	'B'},
		{"noise-cache",	required_argument,	NULL,	'C'},
		{"noise-channel",	required_argument,	NULL,	'D'},
		{"output",		required_argument,	NULL,	'H'},
		{"help",		no_argument,		NULL,	'h'},
		{NULL,	0,	NULL,	0}
	};
	struct options *options = command_line_options_new();

	if(!options)
		return NULL;

	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
	/* snr-cache */
	case 'A':
		options->snr_cache = optarg;
		break;

	/* snr-channel */
	case 'B':
		command_line_set_instrument(options, optarg);
		break;

	/* noise-cache */
	case 'C':
		options->noise_cache = optarg;
		break;

	/* noise-channel */
	case 'D':
		command_line_set_instrument(options, optarg);
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


static COMPLEX16TimeSeries *get_complex16series_from_cache(
	const char *cache_name,
	const char *channel_name
)
{
	char instrument[] = {channel_name[0], channel_name[1], '\0'};
	LALCache *cache;
	LALFrStream *stream;
	COMPLEX8TimeSeries *data;
	COMPLEX16TimeSeries *result;
	unsigned i;

	/* construct stream */
	cache = XLALCacheImport(cache_name);
	if(!cache)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	XLALCacheSieve(cache, 0, 0, instrument, NULL, NULL);
	if(cache->length != 1) {
		XLALPrintError("error: cache must contain exactly 1 file for instrument %s", instrument);
		XLALDestroyCache(cache);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}
	stream = XLALFrStreamCacheOpen(cache);
	XLALDestroyCache(cache);
	if(!stream)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* turn on checking for missing data */
	stream->mode = LAL_FR_STREAM_VERBOSE_MODE;

	/* get data */
	data = XLALFrFileReadCOMPLEX8TimeSeries(stream->file, channel_name, 0);

	/* check for gaps and close */
	XLALFrStreamClose(stream);

	/* error checking */
	if(!data)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* convert to double precision */
	result = XLALCreateCOMPLEX16TimeSeries(data->name, &data->epoch, data->f0, data->deltaT, &data->sampleUnits, data->data->length);
	for(i = 0; i < data->data->length; i++)
		result->data->data[i] = data->data->data[i];
	XLALDestroyCOMPLEX8TimeSeries(data);
	/*for(i = 0; i < result->data->length; i++) fprintf(stderr, "%g+I*%g\n", creal(result->data->data[i]), cimag(result->data->data[i]));*/

	/* done */
	{
	char *s = XLALGPSToStr(NULL, &result->epoch);
	fprintf(stderr, "%s: loaded \"%s\" at %s s, duration=%g s (%d samples)\n", instrument, result->name, s, result->data->length * result->deltaT, result->data->length);
	LALFree(s);
	}
	return result;
}


static COMPLEX16Sequence *get_complex16sequence_from_cache(
	const char *cache_name,
	const char *channel_name
)
{
	char instrument[] = {channel_name[0], channel_name[1], '\0'};
	LALCache *cache;
	LALFrStream *stream;
	COMPLEX16TimeSeries *data;
	COMPLEX16Sequence *result;

	/* construct stream */
	cache = XLALCacheImport(cache_name);
	if(!cache)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	XLALCacheSieve(cache, 0, 0, instrument, NULL, NULL);
	if(cache->length != 1) {
		XLALPrintError("error: cache must contain exactly 1 file for instrument %s", instrument);
		XLALDestroyCache(cache);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}
	stream = XLALFrStreamCacheOpen(cache);
	XLALDestroyCache(cache);
	if(!stream)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* turn on checking for missing data */
	stream->mode = LAL_FR_STREAM_VERBOSE_MODE;

	/* get data */
	data = XLALFrFileReadCOMPLEX16TimeSeries(stream->file, channel_name, 0);

	/* check for gaps and close */
	XLALFrStreamClose(stream);

	/* error checking */
	if(!data)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	{
	char *s = XLALGPSToStr(NULL, &data->epoch);
	fprintf(stderr, "%s: loaded \"%s\" at %s s, duration=%g s (%d samples)\n", instrument, data->name, s, data->data->length * data->deltaT, data->data->length);
	LALFree(s);
	}
	result = data->data;
	data->data = NULL;
	XLALDestroyCOMPLEX16TimeSeries(data);

	/* done */
	return result;
}


/*
 * make the time series span the same intervals
 * FIXME:  this leaves the final GPS times slightly different.  why?
 */


static int time_series_pad(COMPLEX16TimeSeries **series, COMPLEX16Sequence **nseries, int n_series)
{
	int i;
	LIGOTimeGPS start;
	LIGOTimeGPS end;

	start = end = series[0]->epoch;
	XLALGPSAdd(&end, series[0]->deltaT * series[0]->data->length);
	for(i = 1; i < n_series; i++) {
		LIGOTimeGPS this_end = series[i]->epoch;
		XLALGPSAdd(&this_end, series[i]->deltaT * series[i]->data->length);
		if(XLALGPSCmp(&start, &series[i]->epoch) > 0)
			start = series[i]->epoch;
		if(XLALGPSCmp(&end, &this_end) < 0)
			end = this_end;
	}

	for(i = 0; i < n_series; i++) {
		XLALResizeCOMPLEX16TimeSeries(series[i], round(XLALGPSDiff(&start, &series[i]->epoch) / series[i]->deltaT), round(XLALGPSDiff(&end, &start) / series[i]->deltaT));
		XLALResizeCOMPLEX16Sequence(nseries[i], -((int) series[i]->data->length - (int) nseries[i]->length) / 2, series[i]->data->length);
	}

#if 0
	for(i = 0; i < n_series; i++) {
		char *s = XLALGPSToStr(NULL, &series[i]->epoch);
		fprintf(stderr, "zero-padded interval for \"%s\": %s s, duration=%g s (%d samples)\n", series[i]->name, s, series[i]->data->length * series[i]->deltaT, series[i]->data->length);
		LALFree(s);
		{
		unsigned j;
		for(j = 0; j < series[i]->data->length; j++) {
			LIGOTimeGPS t = series[i]->epoch;
			XLALGPSAdd(&t, j * series[i]->deltaT);
			s = XLALGPSToStr(NULL, &t);
			fprintf(stderr, "%s: %g+I*%g\n", s, creal(series[i]->data->data[j]), cimag(series[i]->data->data[j]));
			LALFree(s);
		}
		fprintf(stderr, "\n");
		}
	}
#endif

	return 0;
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
 *                             Projection for antenna
 *
 * ============================================================================
 */


static void FDP(double *fplus, double *fcross, const LALDetector **det, int n, double theta, double phi)
{
	double twopsi;
	double normplus2;
	double normcross2;
	double product;
	int i;

	/* store fp, fc */
	for(i = 0; i < n; i++)
		/* gmst is rotated at the end of this code.
		 * So we can set zero. */
		XLALComputeDetAMResponse(&fplus[i], &fcross[i], det[i]->response, phi, M_PI_2 - theta, 0.0, 0.0);

	/* dominant polarization angle */
	normplus2 = normcross2 = product = 0.0;
	for(i = 0; i < n; i++){
		normplus2 += fplus[i] * fplus[i];
		normcross2 += fcross[i] * fcross[i];
		product += fplus[i] * fcross[i];
	}
	twopsi = atan2(2.0 * product, normplus2 - normcross2) / 2.;

	/* set fp, fc */
	for(i = 0; i < n; i++){
		double temp = cos(twopsi) * fplus[i] + sin(twopsi) * fcross[i];
		fcross[i] = -sin(twopsi) * fplus[i] + cos(twopsi) * fcross[i];
		fplus[i] = temp;
	}

	normplus2 = normcross2 = 0.0;
	for(i = 0; i < n; i++){
		normplus2 += fplus[i] * fplus[i];
		normcross2 += fcross[i] * fcross[i];
	}
	for(i = 0; i < n; i++){
		fplus[i] /= sqrt(normplus2);
		fcross[i] /= sqrt(normcross2);
	}
}


#if 0
static complex double ProjectionMatrix(double theta, double phi, int i, int j, const LALDetector **det, int n)
{
	/* this is general parameterized one */
	double fplus[n], fcross[n];

	FDP(fplus, fcross, det, n, theta, phi);
	return fplus[i] * fplus[j] + fcross[i] * fcross[j];
}
#endif
#if 0
static complex double ProjectionMatrix(double theta, double phi, int i, int j, const LALDetector **det, int n)
{
	/* this is CBC parameterized one for \beta=1 */
	int k;
	double fplus[n], fcross[n];
	double normplus2;
	double normcross2;
	normplus2 = normcross2 = 0.0;

	for(k = 0; k < n; k++){
		/* store fp, fc */
		/* gmst is rotated at the end of this code.
		 * So we can set zero. */
		XLALComputeDetAMResponse(&fplus[k], &fcross[k], det[k]->response, phi, M_PI_2 - theta, 0.0, 0.0);
		/* calculate norms of vector fp & fc */
		normplus2 += fplus[k] * fplus[k];
		normcross2 += fcross[k] * fcross[k];
	}

	/* order between i & j is important. Be consistent with Time delay. */
	return (fplus[i] + I * fcross[i]) * (fplus[j] - I * fcross[j]) / (normplus2 + normcross2);
}
#endif
#if 1
static complex double ProjectionMatrix(double theta, double phi, int i, int j, const LALDetector **det, int n)
{
	/* this is CBC parametrized one for any beta & psi */
	double beta = 1;
	double psi = 0. * M_PI / 6.;

	double fplus[n], fcross[n];
	double normplus2;
	double normcross2;

	int k;
	normplus2 = normcross2 = 0.0;
	for(k = 0; k < n; k++){
		/* store fp, fc */
		/* gmst is rotated at the end of this code.
		 * So we can set zero. */
		XLALComputeDetAMResponse(&fplus[k], &fcross[k], det[k]->response, phi, M_PI_2 - theta, psi, 0.0);
		/* calculate norms of vector fp & fc */
		normplus2 += fplus[k] * fplus[k];
		normcross2 += fcross[k] * fcross[k];
	}

	return (fplus[i] + I * beta * fcross[i]) * (fplus[j] - I * beta * fcross[j]) / (normplus2 + beta*beta * normcross2);
}
#endif


struct ProjectionMatrixWrapperData {
	int i, j;
	const LALDetector **det;
	int n;
};


static complex double ProjectionMatrixWrapper(double theta, double phi, void *_data)
{
	struct ProjectionMatrixWrapperData *data = _data;

	return ProjectionMatrix(theta, phi, data->i, data->j, data->det, data->n);
}


static struct correlator_plan_fd *correlator_plan_mult_by_projection(struct correlator_plan_fd *plan, const struct sh_series *projection)
{
	struct sh_series_array *result;
	struct sh_series_product_plan *product_plan;
	int i;

	result = sh_series_array_new(plan->delay_product->n, plan->delay_product->l_max + projection->l_max, plan->delay_product->series[0].polar && projection->polar);
	if(!result)
		return NULL;
	product_plan = sh_series_product_plan_new(&result->series[0], &plan->delay_product->series[0], projection);
	if(!product_plan) {
		XLALPrintError("sh_series_product_plan_new() failed\n");
		sh_series_array_free(result);
		return NULL;
	}

	for(i = 0; i < plan->delay_product->n; i++) {
		fprintf(stderr, "%.3g%%   \r", 100. * (i + 1.) / plan->delay_product->n);
		if(!sh_series_product(&result->series[i], &plan->delay_product->series[i], projection, product_plan)) {
			XLALPrintError("sh_series_product() failed\n");
			sh_series_array_free(result);
			return NULL;
		}
	}
	fprintf(stderr, "\n");

	sh_series_product_plan_free(product_plan);
	sh_series_array_free(plan->delay_product);
	plan->delay_product = result;

	return plan;
}


static int correlator_network_plan_mult_by_projection(struct correlator_network_plan_fd *plan)
{
	int i;
	struct sh_series *projection = sh_series_new(Projection_lmax, 0);
	const struct instrument_array *instruments = plan->baselines->baselines[0]->instruments;
	const LALDetector **det = malloc(instrument_array_len(instruments) * sizeof(*det));

	if(!projection || !instruments || !det) {
		sh_series_free(projection);
		free(det);
		return -1;
	}

	for(i = 0; i < instrument_array_len(instruments); i++) {
		det[i] = instrument_array_get(instruments, i)->data;
		if(!det[i]) {
			sh_series_free(projection);
			free(det);
			return -1;
		}
	}

	for(i = 0; i < plan->baselines->n_baselines; i++) {
		struct ProjectionMatrixWrapperData data = {
			.i = plan->baselines->baselines[i]->index_a,
			.j = plan->baselines->baselines[i]->index_b,
			.det = det,
			.n = instrument_array_len(instruments)
		};
		if(!sh_series_from_func(projection, ProjectionMatrixWrapper, &data)) {
			sh_series_free(projection);
			free(det);
			return -1;
		}

		/* plot projection */
		char filename[32] ={"\0"};
		char istr[12];
		snprintf(istr, sizeof(istr), "%d", i);
		strcat(filename, "projection");
		strcat(filename, istr);
		strcat(filename, ".fits");
		fprintf(stderr, "generate Projection distribution from %d and %d\n", data.i, data.j);
		if(sh_series_write_healpix_alm(projection, filename)) {
			fprintf(stderr, "write \"%s\" failed\n", filename);
			free(det);
			exit(1);
		}

		/* the delay operator is computed with the baseline rotated
		 * to lie along the z axis to take advantage of the
		 * aziumthal symmetry of that operator, so we need to
		 * rotate the projection function we've just computed the
		 * same way */
		{
		struct sh_series *cpy = sh_series_copy(projection);
		double *R = sh_series_invrot_matrix(plan->baselines->baselines[i]->theta, plan->baselines->baselines[i]->phi);
		struct sh_series_rotation_plan *rot = sh_series_rotation_plan_new(projection, R);
		free(R);
		sh_series_rotate(projection, cpy, rot);
		sh_series_rotation_plan_free(rot);
		sh_series_free(cpy);
		}
		/* multiply each frequency bin of the delay operator by the
		 * projection operator */
		if(!correlator_plan_mult_by_projection(plan->plans[i], projection)) {
			sh_series_free(projection);
			free(det);
			return -1;
		}
		/* replace the correlator's intermediate "1D" storage with
		 * one that does not have azimuthal symmetry, increase it
		 * to the correct harmonic order, and compute a new
		 * rotation plan for it.  NOTE:  the two-stage design that
		 * takes advantage of rotational symmetry about the
		 * baseline axis cannot be used here because after
		 * multiplying by the projection operator above we no
		 * longer have that rotational symmetry;  really we should
		 * compute the delay matrix with the correct orientation
		 * with respect to the sky, and not pay the price of the
		 * rotation */
		if(!sh_series_set_polar(plan->plans[i]->power_1d, 0) || !sh_series_resize(plan->plans[i]->power_1d, plan->plans[i]->delay_product->series[0].l_max)) {
			sh_series_free(projection);
			free(det);
			return -1;
		}
		{
		double *R = sh_series_rot_matrix(plan->baselines->baselines[i]->theta, plan->baselines->baselines[i]->phi);
		if(!R) {
			sh_series_free(projection);
			free(det);
			return -1;
		}
		sh_series_rotation_plan_free(plan->plans[i]->rotation_plan);
		plan->plans[i]->rotation_plan = sh_series_rotation_plan_new(plan->plans[i]->power_1d, R);
		free(R);
		if(!plan->plans[i]->rotation_plan) {
			sh_series_free(projection);
			free(det);
			return -1;
		}
		}
	}

	sh_series_free(projection);
	free(det);

	return 0;
}


/*
 * ============================================================================
 *
 *                                Diagonal part
 *
 * ============================================================================
 */


static int autocorrelator_network_from_projection(struct sh_series *sky, complex double **fseries, struct options *options, unsigned int length)
{
	int i, j;
	struct sh_series *projection = sh_series_new(Projection_lmax, 0);
	const LALDetector **det = malloc(instrument_array_len(options->instruments) * sizeof(*det));

	if(!projection || !det) {
		sh_series_free(projection);
		free(det);
		return -1;
	}

	/* set instruments information */
	struct instrument_array *instruments = instrument_array_new(0);
	for(i = 0; i < instrument_array_len(options->instruments); i++){
		char instrument_name[3] = {options->channels[i][0], options->channels[i][1], '\0'};
		instrument_array_append(instruments, instrument_new_from_name(instrument_name));
	}
	for(i = 0; i < instrument_array_len(instruments); i++) {
		det[i] = instrument_array_get(instruments, i)->data;
		if(!det[i]) {
			sh_series_free(projection);
			free(det);
			return -1;
		}
	}

	for(i = 0; i < instrument_array_len(instruments); i++) {
		/* calc. projection operator */
		struct ProjectionMatrixWrapperData data = {
			.i = i,
			.j = i,
			.det = det,
			.n = instrument_array_len(instruments)
		};
		if(!sh_series_from_func(projection, ProjectionMatrixWrapper, &data)) {
			sh_series_free(projection);
			free(det);
			return -1;
		}

		/* execute calc. */
		double correlator = 0;
		for(j = 0; j < (int) length; j++)
			correlator += fseries[i][j] * conj(fseries[i][j]);
		correlator /= length * length;	// TODO: after considering all TODO, you can decide wheter this line is alive or not.
		for(j = 2; j <= instrument_array_len(instruments); j++)	// TODO: after considering all TODO, you can decide wheter this line is alive or not.
			correlator /= j;
		fprintf(stderr, "diagonal weight %s: %g\n", options->channels[i], correlator);
		sh_series_add(sky, correlator, projection);

#if 1
		/* plot projection */
		char filename[32] ={"\0"};
		char istr[12];
		snprintf(istr, sizeof(istr), "%d", i + instrument_array_len(instruments));
		strcat(filename, "projection");
		strcat(filename, istr);
		strcat(filename, ".fits");
		fprintf(stderr, "generate Projection distribution from %d and %d\n", data.i, data.j);
		sh_series_scale(projection, correlator);
		if(sh_series_write_healpix_alm(projection, filename)) {
			fprintf(stderr, "write \"%s\" failed\n", filename);
			free(det);
			exit(1);
		}
#endif
	}

	sh_series_free(projection);
	return 0;
}


/*
 * ============================================================================
 *
 *                                Whitening
 *
 * ============================================================================
 */


static int threshold(complex double *series, complex double *noise, int length)
{
	int i;

	/* 0.560964 = 10^3 * "numerical error" */
	for(i = 0; i < length; i++)
		if(cabs(*noise++) < 0.560964)
			series[i] = 0;
	/*FIXME: numerical error would depend on emvironment and data length */

	return 0;
}


static int whiten(complex double *series, complex double *noise, int length)
{
	int i;

	for(i = 0; i < length; i++)
		series[i] /= sqrt(cabs(noise[i]));

	return 0;
}

/*
 * ============================================================================
 *
 *                                   prior
 *
 * ============================================================================
 */


static double logprior(double theta, double phi, void *data)
{
	return log(0.25 * sin(theta) / M_PI);
}


static int sh_series_add_logprior(struct sh_series *sky, int lmax)
{
	struct sh_series *sh_series_logprior;
	sh_series_logprior = sh_series_new(lmax, 0);

	if(!sh_series_from_realfunc(sh_series_logprior, logprior, NULL)) {
		sh_series_free(sh_series_logprior);
		return -1;
	}

	sh_series_add(sky, 1.0, sh_series_logprior);
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
	struct options *options;
	COMPLEX16TimeSeries **series;
	COMPLEX16Sequence **nseries;
	struct correlator_network_baselines *baselines;
	complex double **fseries;
	complex double **fnseries;
	fftw_plan *fftplans;
	fftw_plan *nfftplans;
	struct correlator_network_plan_fd *fdplans;
	struct sh_series *sky;
	int k;
	struct timeval t_start, t_end;


	/*
	 * Parse command line.
	 */


	options = command_line_parse(argc, argv);
	if(!options->noise_cache){
		XLALPrintError("need auto-correlation series of consistent template as noise\n");
		exit(1);
	}
	options->instruments->n /= 2;
	/* NOTE: have to free options->channels[k > options->instrument->n].
	 * it's done at the end. */


	/*
	 * Load time series data.
	 */


	series = malloc(instrument_array_len(options->instruments) * sizeof(*series));
	nseries = malloc(instrument_array_len(options->instruments) * sizeof(*nseries));
	fseries = malloc(instrument_array_len(options->instruments) * sizeof(*fseries));
	fnseries = malloc(instrument_array_len(options->instruments) * sizeof(*fnseries));
	fftplans = malloc(instrument_array_len(options->instruments) * sizeof(*fftplans));
	nfftplans = malloc(instrument_array_len(options->instruments) * sizeof(*nfftplans));
	if(!series || !fseries || !fftplans || !nseries|| !fnseries || !nfftplans) {
		XLALPrintError("out of memory\n");
		exit(1);
	}
	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		series[k] = get_complex16series_from_cache(options->snr_cache, options->channels[k]);
		if(!series[k]) {
			XLALPrintError("failure loading snr data\n");
			exit(1);
		}
	}
	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		nseries[k] = get_complex16sequence_from_cache(options->noise_cache, options->channels[k + instrument_array_len(options->instruments)]);	// FIXME: depend a input order from command line, now
		if(!nseries[k]) {
			XLALPrintError("failure loading auto-correlation data\n");
			exit(1);
		}
	}


	/*
	 * Tukey Widow
	 */


	int j;
	for(k = 0; k < instrument_array_len(options->instruments); k++){
		REAL8Window *window = XLALCreateTukeyREAL8Window(series[k]->data->length, 0.1);
		for(j = 0; j < (int) window->data->length; j++){
			series[k]->data->data[j] *= window->data->data[j];
			nseries[k]->data[j] *= window->data->data[j];
		}
	}


	/*
	 * bring to a common interval with zero-padding
	 */


	time_series_pad(series, nseries, instrument_array_len(options->instruments));

#if 0
	// replace with white noise
	for(k = 0; k < (int) series[0]->data->length; k++)
		series[0]->data->data[k] = (double) random() / RAND_MAX + I*(double) random() / RAND_MAX - (0.5 + I*0.5);
	for(k = 1; k < instrument_array_len(options->instruments); k++)
		memcpy(series[k]->data->data, series[0]->data->data, series[0]->data->length * sizeof(*series[0]->data->data));
#endif
	/*for(k = 0; k < instrument_array_len(options->instruments); k++) { unsigned j; for(j = 0; j < series[k]->data->length; j++) fprintf(stderr, "%g+I*%g\n", creal(series[k]->data->data[j]), cimag(series[k]->data->data[j])); fprintf(stderr, "\n"); }*/


	/*
	 * prepare correlator
	 */


	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		complex double *save = malloc(series[k]->data->length * sizeof(*save));
		complex double *nsave = malloc(nseries[k]->length * sizeof(*nsave));
		memcpy(save, series[k]->data->data, series[k]->data->length * sizeof(*save));
		memcpy(nsave, nseries[k]->data, nseries[k]->length * sizeof(*nsave));
		fseries[k] = malloc(series[k]->data->length * sizeof(**fseries));
		fnseries[k] = malloc(nseries[k]->length * sizeof(**fnseries));
		fftplans[k] = correlator_ctseries_to_fseries_plan(series[k]->data->data, fseries[k], series[k]->data->length);
		nfftplans[k] = correlator_ctseries_to_fseries_plan(nseries[k]->data, fnseries[k], nseries[k]->length);
		memcpy(series[k]->data->data, save, series[k]->data->length * sizeof(*save));
		memcpy(nseries[k]->data, nsave, nseries[k]->length * sizeof(*nsave));
		free(save);
		free(nsave);
	}

	fprintf(stderr, "constructing base correlator\n");
	baselines = correlator_network_baselines_new(options->instruments);
	fdplans = correlator_network_plan_fd_new(baselines, series[0]->data->length, series[0]->deltaT);
#if 1
	fprintf(stderr, "applying projection operator\n");
	if(correlator_network_plan_mult_by_projection(fdplans)) {
		fprintf(stderr, "failed\n");
		exit(1);
	}
#endif
	sky = sh_series_new_zero(correlator_network_l_max(baselines, series[0]->deltaT) + Projection_lmax, 0);


	/*
	 * Fourier transform data
	 */


	fprintf(stderr, "starting integration\n");
	gettimeofday(&t_start, NULL);

	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		correlator_ctseries_to_fseries(fftplans[k]);
		correlator_ctseries_to_fseries(nfftplans[k]);
#if 0
		unsigned i; for(i = 0; i < series[k]->data->length; i++) fprintf(stderr, "%g+I*%g\n", creal(fseries[k][i]), cimag(fseries[k][i])); fprintf(stderr, "\n");
		int j = 0;
		double temp = cabs(fseries[k][0]);
		for(i = 1; i < series[k]->data->length; i++) {
			if(temp < cabs(fseries[k][i])){
				j = i;
				temp = cabs(fseries[k][i]);
			}
		}
		fprintf(stderr, "%d %g\n\n", j, temp);
#endif
#if 0
		unsigned i; for(i = 0; i < nseries[k]->length; i++) fprintf(stderr, "%g+I*%g\n", creal(fnseries[k][i]), cimag(fnseries[k][i])); fprintf(stderr, "\n");
		int j = 0;
		double temp = cabs(fnseries[k][0]);
		for(i = 1; i < nseries[k]->length; i++) {
			if(temp < cabs(fnseries[k][i])){
				j = i;
				temp = cabs(fnseries[k][i]);
			}
		}
		fprintf(stderr, "%d %g\n\n", j, temp);
#endif
	}


	/*
	 * Whitening
	 */
	

#if 0
	for(k = 0; k < instrument_array_len(options->instruments); k++){
		threshold(fseries[k], fnseries[k], series[k]->data->length);
		whiten(fseries[k], fnseries[k], series[k]->data->length);
	}
#endif


	/*
	 * Compute angular distribution of integrated cross power.
	 */


	if(!correlator_network_integrate_power_fd(sky, fseries, fdplans)) {
		fprintf(stderr, "cross-correlator failed\n");
		exit(1);
	}
	/* We calculated contrubutions from lower triangular part of Projection
	 * matrix, cross-correlation. However upper one still remain.
	 * Fortunately this calculations are easy because Projection and Time
	 * shift operator is symmetry w.r.t. baseline index, e.g. H1 <--> L1.
	 * Therefore it's OK to twice simply. NOTE: We have to pick up only
	 * real part because correlator is Hermite. However this manipulation
	 * is already done in sh_series_write_healpix_alm(). */
	sh_series_scale(sky, 2.0);
#if 0
	/* add contributions from diagonal part, auto-correlation */
	if(autocorrelator_network_from_projection(sky, fseries, options, series[0]->data->length)){
		fprintf(stderr, "auto-correlator failed\n");
		exit(1);
	}
#endif
	/* multiply baseline numbers. Correlator is normalized by it (See
	 * correlator_network_integrate_power_fd()). However our calculation
	 * doesn't need it. */
	/* multiply data length. Correlator is normalized by it (See
	 * correlator_plan_fd_new() or
	 * correlator_baseline_integrate_power_fd(). Honestly we have either
	 * one enough because Kipp's paper has incorrect calculation. However
	 * codes pass all consistency checks. Threfore our codes doesn't have
	 * error but mistakes. We can neglect the mistakes because our result
	 * is correct). However our Likelihood does't need it. */
	sh_series_scale(sky, fdplans->baselines->n_baselines *4 /** series[0]->data->length*/);


	/*
	 * prior
	 */


#if 1
	fprintf(stderr, "start multiply prior\n");
	if(sh_series_add_logprior(sky, correlator_network_l_max(baselines, series[0]->deltaT) + Projection_lmax)){
		fprintf(stderr, "failur adding prior\n");
		exit(1);
	}
#endif


	/*
	 * Rotate sky.
	 */


	fprintf(stderr, "gmst = %.16g rad\n", gmst_from_epoch_and_offset(series[0]->epoch, series[0]->data->length * series[0]->deltaT / 2.0));
	sh_series_rotate_z(sky, sky, gmst_from_epoch_and_offset(series[0]->epoch, series[0]->data->length * series[0]->deltaT / 2.0));

	gettimeofday(&t_end, NULL);
	fprintf(stderr, "finished integration\n");


	/*
	 * Output.
	 */


	fprintf(stderr, "generate fits file\n");
	if(sh_series_write_healpix_alm(sky, "map.fits")) {
		fprintf(stderr, "write \"map.fits\" failed\n");
		exit(1);
	}


	fprintf(stderr, "analyzed %g s of data in %g s\n", series[0]->data->length * series[0]->deltaT, (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec) * 1e-6);
	fprintf(stderr, "data sample rate was %g Hz\n", 1.0 / series[0]->deltaT);
	fprintf(stderr, "sky was computed to harmonic order l = %d\n", sky->l_max);


	/*
	 * Clean up
	 */


	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		XLALDestroyCOMPLEX16TimeSeries(series[k]);
		XLALDestroyCOMPLEX16Sequence(nseries[k]);
		free(fseries[k]);
		free(fnseries[k]);
		fftw_destroy_plan(fftplans[k]);
		fftw_destroy_plan(nfftplans[k]);
	}
	free(series);
	free(nseries);
	correlator_network_plan_fd_free(fdplans);
	correlator_network_baselines_free(baselines);
	sh_series_free(sky);
	/* below's order must be fixed */
	options->instruments->n *= 2;
	command_line_options_free(options);

	return 0;
}
