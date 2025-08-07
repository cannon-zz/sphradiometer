/*
 * Copyright (C) 2012--2025  Kipp Cannon, Tsutsui Takuya
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
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
#include <lal/LALMalloc.h>
#include <lal/LALSimulation.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
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
	struct instrument_array *instruments;
	char *snr_cache;
	char *noise_cache;
	char **channels;
	char *precalc_path;
	char *psd_cache;
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
		.precalc_path = NULL,
		.psd_cache = NULL,
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
		{"precalc-path",	required_argument,	NULL,	'D'},
		{"psd-cache",	required_argument,	NULL,	'E'},
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

	/* precalc-path */
	case 'D':
		options->precalc_path = optarg;
		break;

	/* psd-cache */
	case 'E':
		options->psd_cache = optarg;
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
	XLALFree(s);
	}
	result = data->data;
	data->data = NULL;
	XLALDestroyCOMPLEX16TimeSeries(data);

	/* done */
	return result;
}


static double *get_PSD_from_cache(
	const char *cache_name,
	const char *channel_name,
	unsigned length
)
{
	char channel[3];
	char instrument[] = {channel_name[0], channel_name[1], '\0'};
	double *data = malloc(length * sizeof(*data));
	FILE *fp = fopen(cache_name, "r");

	if(!fp) {
		fprintf(stderr, "file open error\n");
		return NULL;
	}

	/* count # of characters */
	int n = 0;
	while(fgetc(fp) != EOF)
		n++;
	fseek(fp, 0, SEEK_SET);

	/* read data */
	char filepath[n];
	while(fscanf(fp, "%s %s", channel, filepath) != EOF) {
		if(strcmp(channel, instrument) == 0) {
			unsigned i;
			FILE *fpp = fopen(filepath, "r");
			for(i = 0; i < length; i++)
				fscanf(fpp, "%lf", &data[i]);
			fclose(fpp);
		}
	}

	fclose(fp);
	return data;
}


void PSD_times_sqrt_AutoCorrelation(double *psd, COMPLEX16Sequence *series)
{
	int i;
	complex double *fseries;
	complex double *cpy;
	fftw_plan fftplan;

	/* fourier transform */
	cpy = malloc(series->length * sizeof(*cpy));
	memcpy(cpy, series->data, series->length * sizeof(*cpy));
	fseries = malloc(series->length * sizeof(*fseries));
	fftplan = correlator_ctseries_to_fseries_plan(series->data, fseries, series->length);
	memcpy(series->data, cpy, series->length * sizeof(*cpy));
	free(cpy);
	correlator_ctseries_to_fseries(fftplan);

	/* multiply */
	for(i = 0; i < (int) series->length; i++)
		psd[i] *= sqrt(cabs(fseries[i]));
}


static double *make_one_line(double **mat, int n_det, int length)
{
	int i, j;
	double *vec = malloc(n_det * length * sizeof(*vec));

	for(i = 0; i < n_det; i++)
		for(j = 0; j < length; j++)
			vec[i + n_det * j] = mat[i][j];

	return vec;
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


gsl_matrix_complex *AutoCorrelation2Covariance(COMPLEX16Sequence *nseries, int length)
{
	/* The diagonal components of the covariant matrix must be one.  If the
	 * data length is in even, the center of template auto-correlation
	 * cannot be obtained as 1. */
	if((length & 1) == 0){
		fprintf(stderr, "data length must be odd\n");
		return NULL;
	}

	int i, j;
	gsl_matrix_complex *cov = gsl_matrix_complex_calloc(length, length);

	/* upper half of covariance matrix */
	for(i = 0; i < length / 2; i++)
		for(j = 0; j < (length / 2 + 1 + i); j++)
			gsl_matrix_complex_set(cov, i, j, gsl_complex_rect(creal(nseries->data[length / 2 - i + j]), cimag(nseries->data[length / 2 - i + j])));
	/* lower half of covariance matrix */
	for(i = 0; i < length / 2 + 1; i++)
		for(j = i; j < length; j++)
			gsl_matrix_complex_set(cov, length / 2 + i, j, gsl_complex_rect(creal(nseries->data[- i + j]), cimag(nseries->data[- i + j])));

	return cov;
}


complex double *gsl_matrix_complex_column_series(gsl_matrix_complex *mat, int i)
{
	int j;
	complex double *column = malloc(mat->size1 * sizeof(*column));

	for(j = 0; j < (int) mat->size1; j++)
		column[j] = gsl_matrix_complex_get(mat, j, i);

	return column;
}


double *convert_gsl_vector2array(gsl_vector *vector)
{
	int i;
	double *result = malloc(vector->size * sizeof(*result));

	for(i = 0; i < (int) vector->size; i++)
		result[i] = gsl_vector_get(vector, i);

	return result;
}


complex double **convert_gsl_matrix_complex2columns(gsl_matrix_complex *mat)
{
	/* the index of the input point to the rows, but the index of the
	 * output point to the columns */
	int i;
	complex double **result = malloc(mat->size2 * sizeof(*result));

	for(i = 0; i < (int) mat->size2; i++)
		result[i] = gsl_matrix_complex_column_series(mat, i);

	return result;
}


complex double inner_product(complex double *a, complex double *b, int length)
{
	int i;
	complex double prod = 0;

	for(i = 0; i < length; i++)
		prod += *a++ * conj(*b++);

	return prod;
}


static int add_series(complex double *dest, complex double *a, complex double scale, int length)
{
	int i;

	for(i = 0; i < length; i++)
		dest[i] += a[i] * scale;

	return 0;
}


static int eigens_from_AutoCorrelation(gsl_vector *eval, gsl_matrix_complex *evec, COMPLEX16Sequence *nseries, int length)
{
	/* initialize */
	gsl_matrix_complex *cov = AutoCorrelation2Covariance(nseries, length);
	if(!cov) {
		fprintf(stderr, "failure of making covariant matrix from template auto-correlaiton");
		return 1;
	}
	gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(length);

	/* get eigen systems.  NOTE: The covariance matrix is broken to get
	 * eigenvalues and eigenvectors. */
	gsl_eigen_hermv(cov, eval, evec, w);

	/* free */
	gsl_eigen_hermv_free(w);
	gsl_matrix_complex_free(cov);

	return 0;
}


static int KLwhiten(COMPLEX16TimeSeries *series, COMPLEX16Sequence *nseries)
{
	int i, imax;
	double sumeval;
	double *eval;
	complex double **evec_columns;
	gsl_vector *eval_gsl = gsl_vector_alloc(series->data->length);
	gsl_matrix_complex *evec_gsl = gsl_matrix_complex_alloc(series->data->length, series->data->length);
	complex double *result = calloc(series->data->length, sizeof(*result));

	/* get eigens */
	if(eigens_from_AutoCorrelation(eval_gsl, evec_gsl, nseries, series->data->length)) {
		gsl_matrix_complex_free(evec_gsl);
		gsl_vector_free(eval_gsl);
		return 1;
	}
	gsl_eigen_hermv_sort(eval_gsl, evec_gsl, GSL_EIGEN_SORT_ABS_DESC);

	/* convert gsl to series */
	eval = convert_gsl_vector2array(eval_gsl);
	evec_columns = convert_gsl_matrix_complex2columns(evec_gsl);
	gsl_matrix_complex_free(evec_gsl);
	gsl_vector_free(eval_gsl);

	/* If data does not have numerical errors, then abs(sumeval) ==
	 * series->data->length is satisfied.  Since real data has the
	 * numerical errors, the broken point have to be found. */
	imax = 0;
	sumeval = 0;
	do {
		sumeval += eval[imax++];
	} while(fabs(sumeval) < (double) series->data->length);
	fprintf(stderr, "cutoff index: %d\n", imax - 1);

	/* KL whiten */
	for(i = 0; i < imax; i++) {
		add_series(result, evec_columns[i], inner_product(series->data->data, evec_columns[i], (int) series->data->length) / sqrt(fabs(eval[i])), (int) series->data->length);
	}

	/* store result */
	free(series->data->data);
	series->data->data = result;

	/* free */
	free(eval);
	free(evec_columns);

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
	struct correlator_network_baselines *baselines = NULL;
	struct correlator_network_plan_fd *fdplansp = NULL, *fdplansn = NULL;
	struct autocorrelator_network_plan_fd *fdautoplanp = NULL;
	struct autocorrelator_network_plan_fd *fdautoplann = NULL;
	COMPLEX16TimeSeries **series;
	COMPLEX16Sequence **nseries;
	double **psds = NULL;
	struct sh_series *skyp;
	struct sh_series *skyn;
	struct sh_series *logprior;
	int k;
	struct timeval t_start, t_end;


	/*
	 * Parse command line.
	 */


	options = command_line_parse(argc, argv);


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
		if(!options->noise_cache)
			nseries[k] = convert_TimeSeries2Sequence(get_complex16series_from_cache(options->snr_cache, options->channels[k]));
		else
			nseries[k] = get_complex16sequence_from_cache(options->noise_cache, options->channels[k]);
		if(!nseries[k]) {
			XLALPrintError("failure loading auto-correlation data\n");
			exit(1);
		}
	}


	/*
	 * Whitening
	 */


#if 0
	for(k = 0; k < instrument_array_len(options->instruments); k++)
		if(KLwhiten(series[k], nseries[k]))
			exit(1);
#endif


	preprocess_SNRTimeSeries(series, nseries, instrument_array_len(options->instruments));


	/*
	 * Generate sky alm
	 */


	struct stat statBuf;
	if(instrument_array_len(options->instruments) > 1) {
		if(stat(options->precalc_path, &statBuf)) {
			/* read psds */
			fprintf(stderr, "read psds\n");
			double **temp = malloc(instrument_array_len(options->instruments) * sizeof(*temp));
			if(!temp) {
				XLALPrintError("out of memory\n");
				exit(1);
			}
			for(k = 0; k < instrument_array_len(options->instruments); k++) {
				/* NOTE:  we force the flat PSD code path */
				if(1 || !options->psd_cache) {
					temp[k] = malloc(series[k]->data->length * sizeof(*temp[k]));
					if(!temp[k]) {
						XLALPrintError("failure allocating snr data\n");
						exit(1);
					}
					for(int i = 0; i < (int) series[k]->data->length; i++)
						temp[k][i] = 1.;
				} else {
					temp[k] = get_PSD_from_cache(options->psd_cache, options->channels[k], series[k]->data->length);
					if(!temp[k]) {
						XLALPrintError("failure loading snr data\n");
						exit(1);
					}
				}
				//PSD_times_sqrt_AutoCorrelation(temp[k], nseries[k]);
			}
			psds = transpose_matrix(temp, instrument_array_len(options->instruments), series[0]->data->length);
			for(k = 0; k < instrument_array_len(options->instruments); k++) {
				free(temp[k]);
			}
			free(temp);

			/* prepare pre-calculated objects */
			fprintf(stderr, "constructing base correlator\n");
			baselines = correlator_network_baselines_new(options->instruments);

			logprior = sh_series_log_uniformsky_prior(correlator_network_l_max(baselines, series[0]->deltaT));
			if(!logprior) {
				exit(1);
			}


			fdplansp = correlator_network_plan_fd_new(baselines, series[0]->data->length, series[0]->deltaT);
			fdplansn = correlator_network_plan_fd_copy(fdplansp);
			if(!fdplansp || !fdplansn) {
				fprintf(stderr, "memory error\n");
				exit(1);
			}
#if 1
			fprintf(stderr, "applying projection operator\n");
			if(correlator_network_plan_mult_by_projection(fdplansp, +1, 0, psds)) {
				fprintf(stderr, "positive plan is failed\n");
				exit(1);
			}
			if(correlator_network_plan_mult_by_projection(fdplansn, -1, 0, psds)) {
				fprintf(stderr, "negative plan is failed\n");
				exit(1);
			}
#endif
			fprintf(stderr, "construct plans for auto-correlations\n");
			fdautoplanp = autocorrelator_network_plan_fd_new(fdplansp->baselines->baselines[0]->instruments, +1, 0, psds, (int) fdplansp->plans[0]->delay_product->n, logprior->l_max);
			fdautoplann = autocorrelator_network_plan_fd_new(fdplansn->baselines->baselines[0]->instruments, -1, 0, psds, (int) fdplansn->plans[0]->delay_product->n, logprior->l_max);
			if(!fdautoplanp || !fdautoplann) {
				fprintf(stderr, "memory error\n");
				exit(1);
			}
#if 1
			/* make directories to store pre-calculated objects */
			if(make_precalc_directories(options->precalc_path, instrument_array_len(options->instruments))) {
				fprintf(stderr, "can't make precalculated directories\n");
				exit(1);
			}

			/* save pre-calculated objects */
			fprintf(stderr, "make precalculated objects\n");
			if(write_precalc_time_series_length(series[0]->data->length, options->precalc_path)) {
				fprintf(stderr, "can't save time series length\n");
				exit(1);
			}
			if(write_precalc_logprior(logprior, options->precalc_path)) {
				fprintf(stderr, "false write_precalc_logprior()\n");
				exit(1);
			}
			if(write_precalc_correlator_network_plan_fd(fdplansp, fdplansn, options->precalc_path)) {
				fprintf(stderr, "can't save correlator network plan\n");
				exit(1);
			}
			if(write_precalc_autocorrelator_network_plan_fd(fdautoplanp, fdautoplann, options->precalc_path)) {
				fprintf(stderr, "can't save autocorrelator network plan\n");
				exit(1);
			}
#endif
		} else {
			fprintf(stderr, "read precalculated objects\n");
			unsigned int precalc_length;
			if(read_precalc_time_series_length(&precalc_length, options->precalc_path)) {
				fprintf(stderr, "can't read time series length\n");
				exit(1);
			}
			if(series[0]->data->length != precalc_length) {
				/* NOTE: series[0]->data->length in
				 * precalculated objects must be equivalnt to
				 * series[0]->data->length in this code */
				fprintf(stderr, "time series length is inconsistent with that in precalculated objects.\n");
				exit(1);
			}
			logprior = read_precalc_logprior(options->precalc_path);
			baselines = read_precalc_correlator_network_baselines(options->instruments, options->precalc_path);
			fdplansp = malloc(sizeof(*fdplansp));
			fdplansn = malloc(sizeof(*fdplansn));
			read_precalc_correlator_network_plan_fd(fdplansp, fdplansn, baselines, series[0]->data->length, options->precalc_path);
			fdautoplanp = malloc(sizeof(*fdautoplanp));
			fdautoplann = malloc(sizeof(*fdautoplann));
			read_precalc_autocorrelator_network_plan_fd(fdautoplanp, fdautoplann, options->instruments, series[0]->data->length, options->precalc_path);
		}

		gettimeofday(&t_start, NULL);
		if(generate_alm_skys(&skyp, &skyn, fdplansp, fdplansn, fdautoplanp, fdautoplann, series, nseries, logprior)) {
			fprintf(stderr, "generate_alm_skys error\n");
			exit(1);
		}
	} else {
		if(stat(options->precalc_path, &statBuf)) {
			/* prepare pre-calculated objects */
			logprior = sh_series_log_uniformsky_prior(correlator_baseline_power_l_max_naive(series[0]->data->length * series[0]->deltaT, series[0]->deltaT));
			if(!logprior) {
				exit(1);
			}

#if 1
			/* make directories to store pre-calculated objects */
			if(make_precalc_directories(options->precalc_path, instrument_array_len(options->instruments))) {
				fprintf(stderr, "can't make precalculated directories\n");
				exit(1);
			}

			/* save pre-calculated objects */
			fprintf(stderr, "make precalculated objects\n");
			if(write_precalc_time_series_length(series[0]->data->length, options->precalc_path)) {
				fprintf(stderr, "can't save time series length\n");
				exit(1);
			}
			if(write_precalc_logprior(logprior, options->precalc_path)) {
				fprintf(stderr, "false write_precalc_logprior()\n");
				exit(1);
			}
#endif
		} else {
			fprintf(stderr, "read precalculated objects\n");
			unsigned int precalc_length;
			if(read_precalc_time_series_length(&precalc_length, options->precalc_path)) {
				fprintf(stderr, "can't read time series length\n");
				exit(1);
			}
			if(series[0]->data->length != precalc_length) {
				/* NOTE: series[0]->data->length in
				 * precalculated objects must be equivalnt to
				 * series[0]->data->length in this code */
				fprintf(stderr, "time series length is inconsistent with that in precalculated objects.\n");
				exit(1);
			}
			logprior = read_precalc_logprior(options->precalc_path);
		}

		/* For single detector case, the antenna Projection operator is
		 * unity, that is, independent of the sky direction.  Thus, the
		 * sky probability maps are equal to prior distribution. */
		gettimeofday(&t_start, NULL);
		skyp = sh_series_copy(logprior);
		skyn = sh_series_copy(logprior);
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


	if(psds != NULL) {
		for(k = 0; k < (int) series[0]->data->length; k++) {
			free(psds[k]);
		}
		free(psds);
	}
	for(k = 0; k < instrument_array_len(options->instruments); k++) {
		XLALDestroyCOMPLEX16TimeSeries(series[k]);
		XLALDestroyCOMPLEX16Sequence(nseries[k]);
	}
	free(series);
	free(nseries);
	sh_series_free(skyp);
	sh_series_free(skyn);
	sh_series_free(logprior);
	if(fdplansp != NULL) {
		/* if the # of detectors > 1, then the followings are defined
		 * simultaneously */
		correlator_network_plan_fd_free(fdplansp);
		correlator_network_plan_fd_free(fdplansn);
		autocorrelator_network_plan_fd_free(fdautoplanp);
		autocorrelator_network_plan_fd_free(fdautoplann);
		correlator_network_baselines_free(baselines);
	}
	command_line_options_free(options);

	return 0;
}
