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
#include <sphradiometer/inspiral_sky_map.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/correlator.h>
#include <sphradiometer/diagnostics.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALCache.h>
#include <lal/LALConstants.h>
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
#define DIRECTORY_PATH "/home/tsutsui/Development/sphradiometer"
#define FILE_LEN strlen(DIRECTORY_PATH) + 110


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


static COMPLEX16Sequence *convert_TimeSeries2Sequence(COMPLEX16TimeSeries *series)
{
	return series->data;
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

	long data_length = (long) (round(2 * LAL_REARTH_SI / LAL_C_SI / series[0]->deltaT) + series[0]->data->length);
	for(i = 0; i < n_series; i++) {
		XLALResizeCOMPLEX16TimeSeries(series[i], round(XLALGPSDiff(&start, &series[i]->epoch) / series[i]->deltaT), round(XLALGPSDiff(&end, &start) / series[i]->deltaT));
		XLALResizeCOMPLEX16TimeSeries(series[i], round((series[i]->data->length - data_length) / 2.), data_length);
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
 *                                Whitening
 *
 * ============================================================================
 */


static int threshold(complex double *series, complex double *noise, int length, double threshold)
{
	int i;

	/* 0.560964 = 10^3 * "numerical error" */
	int k = 0;
	for(i = 0; i < length; i++)
		//if(cabs(*noise++) < 0.560964)
		//if(cabs(*noise++) < 1.28791e-13)
		if(cabs(*noise++) < threshold){
			series[i] = 0;
			k += 1;
		}
	/*FIXME: numerical error would depend on emvironment and data length */
	fprintf(stderr, "\nk = %d\nlength = %d\nratio = %g\n\n", k, length, 100. * k / length);

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
 *                                Read precalcs
 *
 * ============================================================================
 */


static struct correlator_baseline *read_precalc_correlator_baseline(const struct instrument_array *instruments, int i)
{
	char filename[FILE_LEN];
	FILE *fp;
	struct correlator_baseline *new = malloc(sizeof(*new));

	/* read no pointer objects */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/baselines/%d/correlator_baseline.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "no correlator_baseline.dat in %d\n", i);
		return NULL;
	}
	/* following order must be consistent with the one in the writer */
	fread(&new->index_a, sizeof(new->index_a), 1, fp);
	fread(&new->index_b, sizeof(new->index_b), 1, fp);
	fread(&new->theta, sizeof(new->theta), 1, fp);
	fread(&new->phi, sizeof(new->phi), 1, fp);
	fclose(fp);

	/* read d */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/baselines/%d/d.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "no d.dat in %d\n", i);
		return NULL;
	}
	new->d = gsl_vector_alloc(3);
	gsl_vector_fread(fp, new->d);
	fclose(fp);

	/* substitute instruments */
	new->instruments = instrument_array_copy(instruments);

	return new;
}


static struct correlator_network_baselines *read_precalc_correlator_network_baselines(const struct instrument_array *instruments)
{
	char filename[FILE_LEN];
	struct correlator_network_baselines *new = malloc(sizeof(*new));
	int i;

	/* read n_baselines */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/n_baselines.dat", DIRECTORY_PATH);
	FILE *fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "no n_baselines.dat\n");
		return NULL;
	}
	fread(&new->n_baselines, sizeof(new->n_baselines), 1, fp);
	fclose(fp);

	/* read baselines */
	new->baselines = malloc(new->n_baselines * sizeof(*new->baselines));
	for(i = 0; i < new->n_baselines; i++)
		new->baselines[i] = read_precalc_correlator_baseline(instruments, i);

	return new;
}


static struct sh_series_rotation_plan *read_precalc_sh_series_rotation_plan(int i)
{
	char filename[FILE_LEN];
	FILE *fp;
	struct sh_series_rotation_plan *new = malloc(sizeof(*new));

	/* read l_max */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/rotation_plan/l_max.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "no %dth l_max.dat\n", i);
		return NULL;
	}
	fread(&new->l_max, sizeof(new->l_max), 1, fp);
	fclose(fp);

	/* read D */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/rotation_plan/D.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "no %dth D.dat\n", i);
		return NULL;
	}
	new->D = malloc(new->l_max * sizeof(*new->D));
	new->D--;
	unsigned int l;
	for(l = 1; l <= new->l_max; l++) {
		new->D[l] = malloc((2 * l + 1) * (2 * l + 1) * sizeof(*new->D[l]));
		fread(new->D[l], sizeof(*new->D[l]), (2 * l + 1) * (2 * l + 1), fp);
		new->D[l] += (2 * l + 1) * l + l;
	}
	fclose(fp);

	return new;
}


static int read_precalc_correlator_plan_fd(struct correlator_plan_fd *fdplanp, struct correlator_plan_fd *fdplann, const struct instrument_array *instruments, int tseries_length, int i)
{
	char filename[FILE_LEN];
	FILE *fp;
	int j;

	if(!fdplanp || !fdplann) {
		fprintf(stderr, "memory of correlator_plan_fd must be allocated\n");
		return -1;
	}

	/* read no pointer objects */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/correlator_plan_fd.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "no %dth correlator_plan_fd.dat\n", i);
		return -1;
	}
	/* following order must be consistent with the one in the writer */
	fread(&fdplanp->delta_t, sizeof(fdplanp->delta_t), 1, fp);
	fread(&fdplanp->transient, sizeof(fdplanp->transient), 1, fp);
	fclose(fp);
	fdplann->delta_t = fdplanp->delta_t;
	fdplann->transient = fdplanp->transient;

	/* read baseline */
	fdplanp->baseline = read_precalc_correlator_baseline(instruments, i);
	fdplann->baseline = correlator_baseline_copy(fdplanp->baseline);

	/* read rotation_plan */
	fdplanp->rotation_plan = read_precalc_sh_series_rotation_plan(i);
	fdplann->rotation_plan = sh_series_rotation_plan_copy(fdplanp->rotation_plan);

	/* read delay_product */
	/* positive case */
	fdplanp->delay_product = malloc(sizeof(*fdplanp->delay_product));
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/delay_product_p/sh_series_array.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "can't open %dth positive sh_series_array.dat\n", i);
		return -1;
	}
	/* following order must be consistent with the one in the writer */
	fread(&fdplanp->delay_product->l_max, sizeof(fdplanp->delay_product->l_max), 1, fp);
	fread(&fdplanp->delay_product->polar, sizeof(fdplanp->delay_product->polar), 1, fp);
	fread(&fdplanp->delay_product->n, sizeof(fdplanp->delay_product->n), 1, fp);
	fread(&fdplanp->delay_product->stride, sizeof(fdplanp->delay_product->stride), 1, fp);
	fdplanp->delay_product->coeff = malloc(fdplanp->delay_product->n * fdplanp->delay_product->stride * sizeof(*fdplanp->delay_product->coeff));
	fread(fdplanp->delay_product->coeff, sizeof(*fdplanp->delay_product->coeff), fdplanp->delay_product->n * fdplanp->delay_product->stride, fp);
	fclose(fp);

	fdplanp->delay_product->series = malloc(fdplanp->delay_product->n * sizeof(*fdplanp->delay_product->series));
	for(j = 0; j < fdplanp->delay_product->n; j++) {
		fdplanp->delay_product->series[j] = (struct sh_series) {
			.l_max = fdplanp->delay_product->l_max,
			.polar = fdplanp->delay_product->polar,
			.coeff = fdplanp->delay_product->coeff + j * fdplanp->delay_product->stride
		};
	}

	/* negative case */
	fdplann->delay_product = malloc(sizeof(*fdplann->delay_product));
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/delay_product_n/sh_series_array.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "rb");
	if(!fp) {
		fprintf(stderr, "can't open %dth negative sh_series_array.dat\n", i);
		return -1;
	}
	/* following order must be consistent with the one in the writer */
	fread(&fdplann->delay_product->l_max, sizeof(fdplann->delay_product->l_max), 1, fp);
	fread(&fdplann->delay_product->polar, sizeof(fdplann->delay_product->polar), 1, fp);
	fread(&fdplann->delay_product->n, sizeof(fdplann->delay_product->n), 1, fp);
	fread(&fdplann->delay_product->stride, sizeof(fdplann->delay_product->stride), 1, fp);
	fdplann->delay_product->coeff = malloc(fdplann->delay_product->n * fdplann->delay_product->stride * sizeof(*fdplann->delay_product->coeff));
	fread(fdplann->delay_product->coeff, sizeof(*fdplann->delay_product->coeff), fdplann->delay_product->n * fdplann->delay_product->stride, fp);
	fclose(fp);

	fdplann->delay_product->series = malloc(fdplann->delay_product->n * sizeof(*fdplann->delay_product->series));
	for(j = 0; j < fdplann->delay_product->n; j++) {
		fdplann->delay_product->series[j] = (struct sh_series) {
			.l_max = fdplann->delay_product->l_max,
			.polar = fdplann->delay_product->polar,
			.coeff = fdplann->delay_product->coeff + j * fdplann->delay_product->stride
		};
	}


	/* read fseries_product */
	/* fseries_product is one of the results.  don't read but initialize. */
	fdplanp->fseries_product = malloc(tseries_length * sizeof(*fdplanp->fseries_product));
	fdplann->fseries_product = malloc(tseries_length * sizeof(*fdplann->fseries_product));

	/* read power_1d */
	/* power_1d is one of the results.  don't read but initialize. */
	fdplanp->power_1d = sh_series_new(fdplanp->rotation_plan->l_max, fdplanp->delay_product->polar);
	fdplann->power_1d = sh_series_new(fdplann->rotation_plan->l_max, fdplann->delay_product->polar);

	return 0;
}


static int read_precalc_correlator_network_plan_fd(struct correlator_network_plan_fd *fdplansp, struct correlator_network_plan_fd *fdplansn, const struct instrument_array *instruments, int tseries_length)
{
	if(!fdplansp || !fdplansn) {
		fprintf(stderr, "memory of correlator_network_plan_fd must be allocated\n");
		return -1;
	}

	/* read baselines */
	fdplansp->baselines = read_precalc_correlator_network_baselines(instruments);
	fdplansn->baselines = correlator_network_baselines_copy(fdplansp->baselines);

	/* read plans */
	fdplansp->plans = malloc(fdplansp->baselines->n_baselines * sizeof(*fdplansp->plans));
	fdplansn->plans = malloc(fdplansn->baselines->n_baselines * sizeof(*fdplansn->plans));
	int i;
	for(i = 0; i < fdplansp->baselines->n_baselines; i++) {
		fdplansp->plans[i] = malloc(sizeof(*fdplansp->plans[i]));
		fdplansn->plans[i] = malloc(sizeof(*fdplansn->plans[i]));
		if(read_precalc_correlator_plan_fd(fdplansp->plans[i], fdplansn->plans[i], instruments, tseries_length, i)) {
			fprintf(stderr, "can't read %dth correlator_plan_fd\n", i);
			return -1;
		}
	}

	return 0;
}


static struct sh_series *read_precalc_logprior(void)
{
	char filename[FILE_LEN];
	sprintf(filename, "%s/precalc/sh_series/logprior.fits", DIRECTORY_PATH);
	return sh_series_read_healpix_alm(filename);
}


/*
 * ============================================================================
 *
 *                                Write precalcs
 *
 * ============================================================================
 */


static int write_precalc_correlator_baseline(const struct correlator_baseline *baseline, int i)
{
	char filename[FILE_LEN];
	FILE *fp;

	/* write no pointer objects */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/baselines/%d/correlator_baseline.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open correlator_baseline.dat\n");
		return -1;
	}
	/* following order must be consistent with the one in the reader */
	fwrite(&baseline->index_a, sizeof(baseline->index_a), 1, fp);
	fwrite(&baseline->index_b, sizeof(baseline->index_b), 1, fp);
	fwrite(&baseline->theta, sizeof(baseline->theta), 1, fp);
	fwrite(&baseline->phi, sizeof(baseline->phi), 1, fp);
	fclose(fp);

	/* write d */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/baselines/%d/d.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open d.dat\n");
		return -1;
	}
	gsl_vector_fwrite(fp, baseline->d);
	fclose(fp);

	/* instruments which is a member of correlator_baseline are given by
	 * command line options.  Then, we don't have to store it. */

	return 0;
}


static int write_precalc_correlator_network_baselines(const struct correlator_network_baselines *baselines)
{
	char filename[FILE_LEN];

	/* write n_baselines */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/n_baselines.dat", DIRECTORY_PATH);
	FILE *fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open n_baselines.dat\n");
		return -1;
	}
	fwrite(&baselines->n_baselines, sizeof(baselines->n_baselines), 1, fp);
	fclose(fp);

	/* write baselines */
	int i;
	for(i = 0; i < baselines->n_baselines; i++) {
		if(write_precalc_correlator_baseline(baselines->baselines[i], i)) {
			fprintf(stderr, "can't save %dth correlator_baseline\n", i);
			return -1;
		}
	}

	return 0;
}


static int write_precalc_sh_series_rotation_plan(const struct sh_series_rotation_plan *plan, int i)
{
	char filename[FILE_LEN];
	FILE *fp;

	/* write l_max */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/rotation_plan/l_max.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open %dth l_max.dat\n", i);
		return -1;
	}
	fwrite(&plan->l_max, sizeof(plan->l_max), 1, fp);
	fclose(fp);

	/* write D */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/rotation_plan/D.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open %dth D.dat\n", i);
		return -1;
	}
	unsigned int l;
	for(l = 1; l <= plan->l_max; l++) {
		plan->D[l] -= (2 * l + 1) * l + l;
		fwrite(plan->D[l], sizeof(*plan->D[l]), (2 * l + 1) * (2 * l + 1), fp);
		plan->D[l] += (2 * l + 1) * l + l;
	}
	fclose(fp);

	return 0;
}


static int write_precalc_correlator_plan_fd(const struct correlator_plan_fd *planp, const struct correlator_plan_fd *plann, int i)
{
	char filename[FILE_LEN];
	FILE *fp;

	/* write no pointer objects */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/correlator_plan_fd.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open %dth correlator_plan_fd.dat\n", i);
		return -1;
	}
	/* following order must be consistent with the one in the reader */
	fwrite(&planp->delta_t, sizeof(planp->delta_t), 1, fp);
	fwrite(&planp->transient, sizeof(planp->transient), 1, fp);
	fclose(fp);

	/* write rotation_plan */
	if(write_precalc_sh_series_rotation_plan(planp->rotation_plan, i)) {
		fprintf(stderr, "can't open %dth rotation_plan.dat\n", i);
		return -1;
	}

	/* write delay_product (sh_series_array) */
	/* positive case */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/delay_product_p/sh_series_array.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open %dth positive sh_series_array.dat\n", i);
		return -1;
	}
	/* following order must be consistent with the one in the reader */
	fwrite(&planp->delay_product->l_max, sizeof(planp->delay_product->l_max), 1, fp);
	fwrite(&planp->delay_product->polar, sizeof(planp->delay_product->polar), 1, fp);
	fwrite(&planp->delay_product->n, sizeof(planp->delay_product->n), 1, fp);
	fwrite(&planp->delay_product->stride, sizeof(planp->delay_product->stride), 1, fp);
	fwrite(planp->delay_product->coeff, sizeof(*planp->delay_product->coeff), planp->delay_product->n * planp->delay_product->stride, fp);
	fclose(fp);

	/* negative case */
	sprintf(filename, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/delay_product_n/sh_series_array.dat", DIRECTORY_PATH, i);
	fp = fopen(filename, "wb");
	if(!fp) {
		fprintf(stderr, "can't open %dth negative sh_series_array.dat\n", i);
		return -1;
	}
	/* following order must be consistent with the one in the reader */
	fwrite(&plann->delay_product->l_max, sizeof(plann->delay_product->l_max), 1, fp);
	fwrite(&plann->delay_product->polar, sizeof(plann->delay_product->polar), 1, fp);
	fwrite(&plann->delay_product->n, sizeof(plann->delay_product->n), 1, fp);
	fwrite(&plann->delay_product->stride, sizeof(plann->delay_product->stride), 1, fp);
	fwrite(plann->delay_product->coeff, sizeof(*plann->delay_product->coeff), plann->delay_product->n * plann->delay_product->stride, fp);
	fclose(fp);

	/* write baseline */
	/* same information as correlator_network_baselines.  Then omit */

	/* write fseries_product */
	/* fseries_product is one of the results.  don't save. */

	/* write power_1d */
	/* power_1d is one of the results.  don't save. */

	return 0;
}


static int write_precalc_correlator_network_plan_fd(const struct correlator_network_plan_fd *fdplansp, const struct correlator_network_plan_fd *fdplansn)
{
	/* write baselines */
	if(write_precalc_correlator_network_baselines(fdplansp->baselines)) {
		fprintf(stderr, "can't save correlator_network_baselines\n");
		return -1;
	}

	/* write plans */
	int i;
	for(i = 0; i < fdplansp->baselines->n_baselines; i++) {
		if(write_precalc_correlator_plan_fd(fdplansp->plans[i], fdplansn->plans[i], i)) {
			fprintf(stderr, "can't save %dth correlator_plan\n", i);
			return -1;
		}
	}

	return 0;
}


static int write_precalc_logprior(const struct sh_series *series)
{
	char filename[FILE_LEN];
	sprintf(filename, "%s/precalc/sh_series/logprior.fits", DIRECTORY_PATH);
	return sh_series_write_healpix_alm(series, filename);
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
	struct correlator_network_baselines *baselines = NULL;
	struct correlator_network_plan_fd *fdplansp, *fdplansn;
	COMPLEX16TimeSeries **series;
	COMPLEX16Sequence **nseries;
	struct sh_series *skyp;
	struct sh_series *skyn;
	struct sh_series *logprior;
	int k;
	struct timeval t_start, t_end;


	/*
	 * Parse command line.
	 */


	options = command_line_parse(argc, argv);
	/* after removing the regulator (avoiding numerical instability), we
	 * need the following check */
	/*
	if(!options->noise_cache){
		XLALPrintError("need auto-correlation series of consistent template as noise\n");
		exit(1);
	}*/
	//options->instruments->n /= 2;	// FIXME: after removing regulator
	/* NOTE: have to free options->channels[k > options->instrument->n].
	 * it's done at the end. */


	/*
	 * Load time series data.
	 */


	series = malloc(instrument_array_len(options->instruments) * sizeof(*series));
	nseries = malloc(instrument_array_len(options->instruments) * sizeof(*nseries));
	if(!series || !nseries) {
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
		//nseries[k] = get_complex16sequence_from_cache(options->noise_cache, options->channels[k + instrument_array_len(options->instruments)]);	// FIXME: depend on a input order from command line, now
		nseries[k] = convert_TimeSeries2Sequence(get_complex16series_from_cache(options->snr_cache, options->channels[k]));	// FIXME: after removing the regularator
		if(!nseries[k]) {
			XLALPrintError("failure loading auto-correlation data\n");
			exit(1);
		}
	}


	/*
	 * Tukey Window
	 */


	for(k = 0; k < instrument_array_len(options->instruments); k++){
		int j;
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
	for(k = 1; k < n; k++)
		memcpy(series[k]->data->data, series[0]->data->data, series[0]->data->length * sizeof(*series[0]->data->data));
#endif
	/*for(k = 0; k < n; k++) { unsigned j; for(j = 0; j < series[k]->data->length; j++) fprintf(stderr, "%g+I*%g\n", creal(series[k]->data->data[j]), cimag(series[k]->data->data[j])); fprintf(stderr, "\n"); }*/
#if 0
	/* replace with exp(2*pi*i*x) for determination of numerical error */
	for(k = 0; k < n; k++) {
		unsigned j;
		for(j = 0; j < series[k]->data->length; j++)
			series[k]->data->data[j] = cos(2 * M_PI * j / series[k]->data->length) + I * sin(2 * M_PI * j / series[k]->data->length);
	}
#endif


	/*
	 * Generate sky alm
	 */


	char dirname[FILE_LEN];
	struct stat statBuf;
	sprintf(dirname, "%s/precalc", DIRECTORY_PATH);
	if(stat(dirname, &statBuf)) {
		/* make directories to store pre-calculated objects */
		mode_t mode = S_IRWXU | S_IRWXG | S_IRWXO;

		fprintf(stderr, "make precalc directory in %s\n", dirname);
		if(mkdir(dirname, mode)) {
			fprintf(stderr, "error making %s\n", dirname);
			exit(1);
		}
		sprintf(dirname, "%s/precalc/correlator_network_plan_fd", DIRECTORY_PATH);
		if(mkdir(dirname, mode)) {
			fprintf(stderr, "error making %s\n", dirname);
			exit(1);
		}
		sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines", DIRECTORY_PATH);
		if(mkdir(dirname, mode)) {
			fprintf(stderr, "error making %s\n", dirname);
			exit(1);
		}
		sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/baselines", DIRECTORY_PATH);
		if(mkdir(dirname, mode)) {
			fprintf(stderr, "error making %s\n", dirname);
			exit(1);
		}
		for(k = 0; k < instrument_array_len(options->instruments); k++){
			sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_network_baselines/baselines/%d", DIRECTORY_PATH, k);
			if(mkdir(dirname, mode)) {
				fprintf(stderr, "error making %s\n", dirname);
				exit(1);
			}
		}
		sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd", DIRECTORY_PATH);
		if(mkdir(dirname, mode)) {
			fprintf(stderr, "error making %s\n", dirname);
			exit(1);
		}
		for(k = 0; k < instrument_array_len(options->instruments); k++){
			sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d", DIRECTORY_PATH, k);
			if(mkdir(dirname, mode)) {
				fprintf(stderr, "error making %s\n", dirname);
				exit(1);
			}
			sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/delay_product_n", DIRECTORY_PATH, k);
			if(mkdir(dirname, mode)) {
				fprintf(stderr, "error making %s\n", dirname);
				exit(1);
			}
			sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/delay_product_p", DIRECTORY_PATH, k);
			if(mkdir(dirname, mode)) {
				fprintf(stderr, "error making %s\n", dirname);
				exit(1);
			}
			sprintf(dirname, "%s/precalc/correlator_network_plan_fd/correlator_plan_fd/%d/rotation_plan", DIRECTORY_PATH, k);
			if(mkdir(dirname, mode)) {
				fprintf(stderr, "error making %s\n", dirname);
				exit(1);
			}
		}
		sprintf(dirname, "%s/precalc/sh_series", DIRECTORY_PATH);
		if(mkdir(dirname, mode)) {
			fprintf(stderr, "error making %s\n", dirname);
			exit(1);
		}

		/* prepare pre-calculated objects */
		fprintf(stderr, "constructing base correlator\n");
		baselines = correlator_network_baselines_new(options->instruments);

		logprior = sh_series_log_uniformsky_prior(correlator_network_l_max(baselines, series[0]->deltaT) + Projection_lmax);
		if(!logprior) {
			exit(1);
		}


		fdplansp = correlator_network_plan_fd_new(baselines, series[0]->data->length, series[0]->deltaT);
		fdplansn = correlator_network_plan_fd_copy(fdplansp);
#if 1
		fprintf(stderr, "applying projection operator\n");
		if(correlator_network_plan_mult_by_projection(fdplansp, +1, 0)) {
			fprintf(stderr, "positive plan is failed\n");
			exit(1);
		}
		if(correlator_network_plan_mult_by_projection(fdplansn, -1, 0)) {
			fprintf(stderr, "negative plan is failed\n");
			exit(1);
		}
#endif
#if 1
		/* save pre-calculated objects */
		fprintf(stderr, "make precalculated objects\n");
		if(write_precalc_logprior(logprior))
			fprintf(stderr, "false write_precalc_logprior()\n");
		if(write_precalc_correlator_network_plan_fd(fdplansp, fdplansn)) {
			fprintf(stderr, "can't save positive plan\n");
			exit(1);
		}
#endif
	} else {
		/* NOTE: series[0]->data->length in precalculated objects must be
		 * equivalnt to series[0]->data->length in this code */
		fprintf(stderr, "read precalculated objects\n");
		logprior = read_precalc_logprior();
		fdplansp = malloc(sizeof(*fdplansp));
		fdplansn = malloc(sizeof(*fdplansn));
		read_precalc_correlator_network_plan_fd(fdplansp, fdplansn, options->instruments, series[0]->data->length);
	}


	gettimeofday(&t_start, NULL);
	if(generate_alm_skys(&skyp, &skyn, fdplansp, fdplansn, series, nseries, logprior)) {
		fprintf(stderr, "generate_alm_skys error\n");
		exit(1);
	}


	/*
	 * Output.
	 */


	fprintf(stderr, "generate fits file\n");
	if(sh_series_write_healpix_alm(skyp, "coeffp.fits")) {
		fprintf(stderr, "write \"coeffp.fits\" failed\n");
		exit(1);
	}
	if(sh_series_write_healpix_alm(skyn, "coeffn.fits")) {
		fprintf(stderr, "write \"coeffn.fits\" failed\n");
		exit(1);
	}


	gettimeofday(&t_end, NULL);
	fprintf(stderr, "analyzed %g s of data in %g s\n", series[0]->data->length * series[0]->deltaT, (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec) * 1e-6);
	fprintf(stderr, "data sample rate was %g Hz\n", 1.0 / series[0]->deltaT);
	fprintf(stderr, "sky was computed to harmonic order l = %d\n", skyp->l_max);


	/*
	 * Clean up
	 */


	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		XLALDestroyCOMPLEX16TimeSeries(series[k]);
		XLALDestroyCOMPLEX16Sequence(nseries[k]);
	}
	free(series);
	free(nseries);
	sh_series_free(skyp);
	sh_series_free(skyn);
	sh_series_free(logprior);
	correlator_network_plan_fd_free(fdplansp);
	correlator_network_plan_fd_free(fdplansn);
	if(baselines != NULL)
		correlator_network_baselines_free(baselines);
	/* following order must be fixed.
	 * all free except for option must be above */
	//options->instruments->n *= 2;	// FIXME: after removing regulator
	command_line_options_free(options);

	return 0;
}
