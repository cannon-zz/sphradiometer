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
#include <getopt.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/correlator.h>
#include <diagnostics.h>
#include <instruments.h>
#include <chealpix.h>
#include <hdf5.h>

#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimulation.h>
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
	char *data1_cache_name;
	char *data1_channel_name;
	char *data2_cache_name;
	char *data2_channel_name;
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
		.instruments = instrument_array_new(2),
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

	instrument_array_set(options->instruments, n, instrument_from_LALDetector(XLALDetectorPrefixToLALDetector(instrument_name)));

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
		{"data1-cache-name",	required_argument,	NULL,	'A'},
		{"data1-channel-name",	required_argument,	NULL,	'B'},
		{"data2-cache-name",	required_argument,	NULL,	'a'},
		{"data2-channel-name",	required_argument,	NULL,	'b'},
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


static COMPLEX8TimeSeries *get_complex8series_from_cache(
	const char *cache_name,
	const char *channel_name,
	LIGOTimeGPS start,
	double duration
)
{
	LALCache *cache;
	LALFrStream *stream;
	COMPLEX8TimeSeries *data;
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
	data = XLALFrStreamReadCOMPLEX8TimeSeries(stream, channel_name, &start, duration, 0);

	/* check for gaps and close */
	gap = stream->state & LAL_FR_STREAM_GAP;
	XLALFrStreamClose(stream);

	/* error checking */
	if(!data)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if(gap) {
		XLALDestroyCOMPLEX8TimeSeries(data);
		XLALPrintError("error: gap detected in input data");
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/* done */
	return data;
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


/* if normalization is needed, set polarization = 1.
 * if it is NOT needed, set polarization = 0. */
static void FDP(double *fplus, double *fcross, const LALDetector **det, int n, double theta, double phi, int normalization)
{
	double twopsi;
	double normplus2;
	double normcross2;
	double product;
	int i;

	// store fp, fc
	for(i = 0; i < n; i++)
		XLALComputeDetAMResponse(&fplus[i], &fcross[i], det[i]->response, phi, M_PI_2 - theta, 0.0, 0.0);

	// dominant polarization angle
	normplus2 = normcross2 = product = 0.0;
	for(i = 0; i < n; i++){
		normplus2 += fplus[i] * fplus[i];
		normcross2 += fcross[i] * fcross[i];
		product += fplus[i] * fcross[i];
	}
	twopsi = atan(2.0 * product / (normplus2 - normcross2)) / 2.;

	// normalization if necessary
	if(normalization)
		for(i = 0; i < n; i++){
			fplus[i] /= sqrt(normplus2);
			fcross[i] /= sqrt(normcross2);
		}

	// set fp, fc
	for(i = 0; i < n; i++){
		double temp = cos(twopsi) * fplus[i] + sin(twopsi) * fcross[i];
		fcross[i] = -sin(twopsi) * fplus[i] + cos(twopsi) * fcross[i];
		fplus[i] = temp;
	}
}


static double ProjectionMatrix(double theta, double phi, int i, int j, const LALDetector **det, int n)
{
	double fplus[n], fcross[n];

	FDP(fplus, fcross, det, n, theta, phi, 1);
	return fplus[i] * fplus[j] + fcross[i] * fcross[j];
}


struct ProjectionMatrixWrapperData {
	int i, j;
	const LALDetector **det;
	int n;
};

static double ProjectionMatrixWrapper(double theta, double phi, void *_data)
{
	struct ProjectionMatrixWrapperData *data = _data;

	return ProjectionMatrix(theta, phi, data->i, data->j, data->det, data->n);
}


static struct correlator_plan_fd *correlator_plan_mult_by_projection(struct correlator_plan_fd *plan, const struct sh_series *projection)
{
	struct sh_series *result;
	struct sh_series_product_plan *product_plan;
	int i;

	result = sh_series_new(plan->delay_product->l_max, 0);
	if(!result)
		return NULL;
	product_plan = sh_series_product_plan_new(result, &plan->delay_product->series[0], projection);
	if(!product_plan) {
		sh_series_free(result);
		return NULL;
	}

	for(i = 0; i < plan->delay_product->n; i++) {
		sh_series_product(result, &plan->delay_product->series[i], projection, product_plan);
		sh_series_assign(&plan->delay_product->series[i], result);
	}

	sh_series_product_plan_free(product_plan);
	sh_series_free(result);

	return plan;
}


static struct correlator_network_plan_fd *correlator_network_plan_mult_by_projection(struct correlator_network_plan_fd *plan)
{
	int i;
	const struct instrument_array *instruments = plan->baselines->baselines[0]->instruments;
	const LALDetector **det = malloc(instrument_array_len(instruments) * sizeof(*det));

	for(i = 0; i < instrument_array_len(instruments); i++)
		det[i] = instrument_array_get(instruments, i)->data;

	for(i = 0; i < plan->baselines->n_baselines; i++) {
		struct sh_series *projection = sh_series_new(8, 0);
		struct ProjectionMatrixWrapperData data = {
			.i = plan->baselines->baselines[i]->index_a,
			.j = plan->baselines->baselines[i]->index_b,
			.det = det,
			.n = instrument_array_len(instruments)
		};
		sh_series_from_realfunc(projection, ProjectionMatrixWrapper, &data);
		correlator_plan_mult_by_projection(plan->plans[i], projection);
		sh_series_free(projection);
	}

	free(det);
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
	//struct options *options;
	COMPLEX8TimeSeries *series[2];
	struct correlator_network_baselines *baselines;
	double *tseries;
	REAL8Window *window;
	complex double *fseries[2];
	fftw_plan fftplans[2];
	struct correlator_network_plan_fd *fdplans;
	int sky_l_max;
	struct sh_series *sky;
	int k;
	struct timeval t_start, t_end;
	int start_sample, time_series_length;
	double gmst;
	const LALDetector det[2];
	int i;
	hid_t file_id, dataset_id, dataspace_id, filespace_id; //
	herr_t status;
	double *dset_data;
	int rank;
	hsize_t dims[2];


	// FIXME: set tseries and initial_time
//	file_id = H5Fopen("output.h5", H5F_ACC_RDWR, H5P_DEFAULT);	// open an existing file
//	dataset_id = H5Dopen(file_id, "/dset", H5P_DEFAULT);		// open an existing dataset
//	filespace_id = H5Dget_space(dataset_id);			// read the dataset
//	rank = H5Sget_simple_extent_ndims(filespace_id);		//
//	status = H5Sget_simple_extent_dims(filespace_id, dims, NULL);	//
//	dset_data = calloc(sizeof(double), dims[0]*dims[1]);		//
//	status = H5Dread(dataset_id, H5T_NAIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data); //
	

	complex float tseries_h1[] = {-3.983366012573242 -I*1.126609563827515 ,
		-4.169796466827393 -I*0.1484837979078293 ,
		-3.969766616821289+ I*0.7086409330368042 ,
		-3.590789794921875+ I*1.213011264801025 ,
		-3.466731309890747+ I*1.413019061088562 ,
		-3.683328866958618+ I*1.755069017410278 ,
		-3.960615873336792+ I*2.538711786270142 ,
		-3.926766633987427+ I*3.765216112136841 ,
		-3.297438383102417+ I*5.179051399230957 ,
		-1.979390144348145+ I*6.397994041442871 ,
		-0.1174830719828606+ I*7.072439193725586 ,
		1.977809309959412+ I*6.977499485015869 ,
		3.828668355941772+ I*6.035830974578857 ,
		5.097399234771729+ I*4.576348781585693 ,
		5.724094867706299+ I*2.90929913520813 ,
		5.729362964630127+ I*1.364532113075256 ,
		5.399999618530273+ I*0.1394371688365936 ,
		4.840004920959473 -I*0.8327525854110718 ,
		4.231722354888916 -I*1.482487440109253 ,
		3.69499659538269 -I*1.961324572563171 ,
		3.169871807098389 -I*2.357835292816162 ,
	};
	complex float tseries_l1[] = {4.783010959625244+ I*4.272365093231201 ,
		5.034692287445068+ I*3.548733711242676 ,
		5.366576194763184+ I*3.200068950653076 ,
		5.998232364654541+ I*3.035331010818481 ,
		7.141907691955566+ I*2.736834764480591 ,
		8.686773300170898+ I*1.842572093009949 ,
		10.26880168914795+ I*0.006665350869297981 ,
		11.27248287200928 -I*2.897961139678955 ,
		11.00472545623779 -I*6.521040439605713 ,
		9.039669036865234 -I*10.08332538604736 ,
		5.528451919555664 -I*12.57289409637451 ,
		1.317213177680969 -I*13.27334213256836 ,
		-2.55349063873291 -I*12.10455894470215 ,
		-5.085619449615479 -I*9.672429084777832 ,
		-6.066431045532227 -I*7.106531620025635 ,
		-6.081847667694092 -I*5.159073829650879 ,
		-5.780312061309814 -I*4.009140491485596 ,
		-5.821097850799561 -I*3.332641124725342 ,
		-6.196267127990723 -I*2.49298095703125 ,
		-6.42664098739624 -I*1.259080767631531 ,
		-6.144642353057861+ I*0.1295482814311981
	};
	complex float tseries_v1[] = {-0.7310531735420227+ I*1.292859077453613 ,
		-0.5984863638877869+ I*1.589023351669312 ,
		-0.3703035116195679+ I*1.957797408103943 ,
		0.06488069146871567+ I*2.389910221099854 ,
		0.8175532817840576+ I*2.72331714630127 ,
		1.821033835411072+ I*2.719309568405151 ,
		2.848188877105713+ I*2.27576470375061 ,
		3.702071189880371+ I*1.421364903450012 ,
		4.21214771270752+ I*0.230946958065033 ,
		4.220407485961914 -I*1.115646958351135 ,
		3.713864088058472 -I*2.333185911178589 ,
		2.857001543045044 -I*3.230406999588013 ,
		1.830233216285706 -I*3.75080394744873 ,
		0.7930911183357239 -I*3.883560180664062 ,
		-0.08061070740222931 -I*3.709914207458496 ,
		-0.7251806855201721 -I*3.417925357818604 ,
		-1.229905009269714 -I*3.130890607833862 ,
		-1.694764614105225 -I*2.844826221466064 ,
		-2.152580738067627 -I*2.520911693572998 ,
		-2.606910943984985 -I*2.119481086730957 ,
		-3.038182497024536 -I*1.595876097679138
	};


	/*
	 * Parse command line.
	 */


	//options = command_line_parse(argc, argv);
	struct options *options = command_line_options_new();
	options->integration_duration = 21 * 0.00048828125;
	options->instruments->instruments[0] = instrument_from_LALDetector(XLALDetectorPrefixToLALDetector("H1"));
	options->instruments->instruments[1] = instrument_from_LALDetector(XLALDetectorPrefixToLALDetector("L1"));
	//options->instruments[2] = instrument_from_LALDetector(XLALDetectorPrefixToLALDetector("V1"));
	time_series_length = 21; //


	/*
	 * Load time series data.
	 */


	/*
	series[0] = get_complex8series_from_cache(options->data1_cache_name, options->data1_channel_name, options->analysis_start, options->analysis_duration);
	series[1] = get_complex8series_from_cache(options->data2_cache_name, options->data2_channel_name, options->analysis_start, options->analysis_duration);

	if(!series[0] || !series[1]) {
		XLALPrintError("failure loading data\n");
		exit(1);
	}
	if(series[0]->deltaT != series[1]->deltaT) {
		fprintf(stderr, "sample rate mismatch\n");
		exit(1);
	}
	*/

	LIGOTimeGPS epoch1, epoch2;
	XLALINT8NSToGPS(&epoch1, 0);
	XLALINT8NSToGPS(&epoch2, 1);
	series[0] = XLALCreateCOMPLEX8TimeSeries("H1", &epoch1, 0, 0.00048828125, &lalSecondUnit, 21);
	free(series[0]->data->data);
	series[0]->data->data = tseries_h1;

	series[1] = XLALCreateCOMPLEX8TimeSeries("L1", &epoch2, 0, 0.00048828125, &lalSecondUnit, 21);
	free(series[1]->data->data);
	series[1]->data->data = tseries_l1;

	window = XLALCreateRectangularREAL8Window(options->integration_duration / series[0]->deltaT);
	if(window->data->length * series[0]->deltaT != options->integration_duration) {
		fprintf(stderr, "integration duration not an integer number of samples\n");
		exit(1);
	}

	tseries = malloc(window->data->length * sizeof(*tseries));
	fseries[0] = malloc(window->data->length * sizeof(**fseries));
	fseries[1] = malloc(window->data->length * sizeof(**fseries));

	fftplans[0] = correlator_tseries_to_fseries_plan(tseries, fseries[0], window->data->length);
	fftplans[1] = correlator_tseries_to_fseries_plan(tseries, fseries[1], window->data->length);


	/*
	 * prepare correlator
	 */


	baselines = correlator_network_baselines_new(options->instruments);
	fdplans = correlator_network_plan_fd_new(baselines, window->data->length, series[0]->deltaT);
	sky_l_max = correlator_network_l_max(baselines, series[0]->deltaT);
	sky = sh_series_new_zero(sky_l_max, 0);


	/*
	 * Loop over integration chunk.
	 */


	fprintf(stderr, "%g", time_series_length);
	fprintf(stderr, "starting integration\n");
	gettimeofday(&t_start, NULL);
	for(start_sample = 0; start_sample + time_series_length <= series[0]->data->length; start_sample += time_series_length / 2) {
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

	status = H5Dclose(dataset_id);	// close the dataset
	status = H5Fclose(file_id);	// close the file
	free(dset_data);		//

	return 0;
}
