/*
 * Copyright (C) 2020  Takuya Tsutsui
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
#include <math.h>
#include <sphradiometer/inspiral_sky_map.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/correlator.h>
#include <stdlib.h>
#include <string.h>

#include <lal/Date.h>
#include <lal/DetResponse.h>
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
		/* gmst is rotated in generate_alm_sky().
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


static complex double ExcessProjectionMatrix(double theta, double phi, int i, int j, const LALDetector **det, int n)
{
	/* this is general parameterized one */
	double fplus[n], fcross[n];

	FDP(fplus, fcross, det, n, theta, phi);
	return fplus[i] * fplus[j] + fcross[i] * fcross[j] - (int)(i == j);
}


static complex double CBCProjectionMatrix(double theta, double phi, int i, int j, const LALDetector **det, int n, double beta, double psi)
{
	/* this is CBC parameterized one for arbitrary beta & psi */
	double fplus[n], fcross[n];
	double normplus2;
	double normcross2;

	int k;
	normplus2 = normcross2 = 0.0;
	for(k = 0; k < n; k++){
		/* store fp, fc */
		/* gmst is rotated in generate_alm_sky().
		 * So we can set zero. */
		XLALComputeDetAMResponse(&fplus[k], &fcross[k], det[k]->response, phi, M_PI_2 - theta, psi, 0.0);
		/* calculate norms of vector fp & fc */
		normplus2 += fplus[k] * fplus[k];
		normcross2 += fcross[k] * fcross[k];
	}

	return (fplus[i] + I * beta * fcross[i]) * (fplus[j] - I * beta * fcross[j]) / (normplus2 + beta*beta * normcross2) - (int)(i == j);
}


struct ProjectionMatrixWrapperData {
	int i, j;
	const LALDetector **det;
	int n;
	double beta;
	double psi;
};


static complex double ProjectionMatrixWrapper(double theta, double phi, void *_data)
{
	struct ProjectionMatrixWrapperData *data = _data;

#if 0
	return ExcessProjectionMatrix(theta, phi, data->i, data->j, data->det, data->n);
#else
	return CBCProjectionMatrix(theta, phi, data->i, data->j, data->det, data->n, data->beta, data->psi);
#endif
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


int correlator_network_plan_mult_by_projection(struct correlator_network_plan_fd *plan, double beta, double psi)
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
			.n = instrument_array_len(instruments),
			.beta = beta,
			.psi = psi
		};
		if(!sh_series_from_func(projection, ProjectionMatrixWrapper, &data)) {
			sh_series_free(projection);
			free(det);
			return -1;
		}

#if 0
		/* plot projection */
		char filename[32] ={"\0"};
		char istr[12];
		snprintf(istr, sizeof(istr), "%d", i);
		strcat(filename, "projection");
		strcat(filename, istr);
		strcat(filename, ".fits");
		if(sh_series_write_healpix_alm(projection, filename)) {
			fprintf(stderr, "write \"%s\" failed\n", filename);
			free(det);
			exit(1);
		}
#endif
		fprintf(stderr, "generate Projection distribution from %d and %d\n", data.i, data.j);

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
 *                                Diagonal Part
 *
 * ============================================================================
 */


static int autocorrelator_network_from_projection(struct sh_series *sky, complex double **fseries, struct instrument_array *instruments, unsigned int length)
{
	// FIXME: This is broken
	int i, j;
	struct sh_series *projection = sh_series_new(Projection_lmax, 0);
	const LALDetector **det = malloc(instrument_array_len(instruments) * sizeof(*det));

	if(!projection || !det) {
		sh_series_free(projection);
		free(det);
		return -1;
	}

	/* set instruments information */
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
		correlator /= length * length;	// TODO: after considering all TODO, you can decide whether this line is alive or not.
		for(j = 2; j <= instrument_array_len(instruments); j++)
			correlator /= j;
		fprintf(stderr, "diagonal weight for channel %d: %g\n", i, correlator);
		sh_series_add(sky, correlator, projection);

#if 0
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


/*
 * ============================================================================
 *
 *                                   prior
 *
 * ============================================================================
 */


static double log_uniformsky_prior(double theta, double phi, void *data)
{
	return log(0.25 * sin(theta) / M_PI);
}


struct sh_series *sh_series_log_uniformsky_prior(int l_max)
{
	struct sh_series *new = sh_series_new(l_max, 0);
	if(!new)
		return NULL;
	if(!sh_series_from_realfunc(new, log_uniformsky_prior, NULL)) {
		sh_series_free(new);
		return NULL;
	}
	return new;
}


/*
 * ============================================================================
 *
 *                                  alm map
 *
 * ============================================================================
 */



static int generate_alm_sky(struct sh_series *sky, struct correlator_network_plan_fd *fdplans, complex double **fseries, complex double **fnseries, double gmst, struct sh_series *logprior)
{
	/*
	 * Compute angular distribution of integrated cross power.
	 */


	fprintf(stderr, "starting integration\n");
	if(!correlator_network_integrate_power_fd(sky, fseries, fdplans)) {
		fprintf(stderr, "cross-correlator failed\n");
		return -1;
	}

	/* We calculated contrubutions from lower triangular part of Projection
	 * matrix, cross-correlation. However upper one still remain.
	 * Fortunately this calculations are easy because Projection and Time
	 * shift operator is symmetric w.r.t. baseline index, e.g. H1 <--> L1.
	 * Therefore it's OK to twice simply. NOTE: We have to pick up only
	 * real part because correlator is Hermite. However this manipulation
	 * is already done in sh_series_write_healpix_alm(). */
	//sh_series_scale(sky, 2.0);
	/* above one is for only cross-correlator, below one for cross- and
	 * auto-correlator */
#if 0
	/* add contributions from diagonal part, auto-correlation. FIXME: check
	 * consistency between above and below sh_series_scale() */
	if(autocorrelator_network_from_projection(sky, fseries, instruments, series[0]->data->length)){
		fprintf(stderr, "auto-correlator failed\n");
		exit(1);
	}
#endif
	/* multiply baseline numbers. Correlator is normalized by it (See
	 * correlator_network_integrate_power_fd()). However our calculation
	 * doesn't need it. */
	/* Correlator is normalized by data length (See
	 * correlator_plan_fd_new() or
	 * correlator_baseline_integrate_power_fd(). Honestly we have either
	 * one enough because Kipp's paper has incorrect calculation. However
	 * codes pass all consistency checks. Threfore our codes doesn't have
	 * error but mistakes. We can neglect the mistakes because our result
	 * is correct). However our Likelihood does't need it. */
	/* Now our calculation has no whitening. To justify it, we recognize it
	 * as multiplying a regulator. The regulator should be normalized. One
	 * SNR timeseries give us one 2 as normalized factor. Hence multiply 4.
	 * */
	sh_series_scale(sky, fdplans->baselines->n_baselines * 8.0);
	fprintf(stderr, "finished integration\n");


	/*
	 * prior
	 */


#if 1
	fprintf(stderr, "start multiply prior\n");
	sh_series_add(sky, 1.0, logprior);
#endif


	/*
	 * Rotate sky.
	 */


	fprintf(stderr, "gmst = %.16g rad\n", gmst);
	sh_series_rotate_z(sky, sky, gmst);

	return 0;
}


int generate_alm_skys(struct sh_series **skyp, struct sh_series **skyn, struct correlator_network_plan_fd *fdplansp, struct correlator_network_plan_fd *fdplansn, COMPLEX16TimeSeries **series, COMPLEX16Sequence **nseries, struct sh_series *logprior)
{
	complex double **fseries;
	complex double **fnseries;
	fftw_plan *fftplans;
	fftw_plan *nfftplans;
	int n = instrument_array_len(fdplansp->baselines->baselines[0]->instruments);
	int k;


	/*
	 * Prepare some serieses
	 */


	fseries = malloc(n * sizeof(*fseries));
	fnseries = malloc(n * sizeof(*fnseries));
	fftplans = malloc(n * sizeof(*fftplans));
	nfftplans = malloc(n * sizeof(*nfftplans));
	if(!fseries || !fftplans || !fnseries || !nfftplans) {
		XLALPrintError("out of memory\n");
		return -1;
	}


	/*
	 * prepare correlator
	 */


	for(k = 0; k < n; k++) {
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


	/*
	 * Fourier transform data
	 */


	for(k = 0; k < n; k++) {
		correlator_ctseries_to_fseries(fftplans[k]);
		correlator_ctseries_to_fseries(nfftplans[k]);
#if 0
		unsigned i; for(i = 0; i < series[k]->data->length; i++) fprintf(stderr, "%3d  %g  %g\n", i, cabs(fseries[k][i]), cabs(fnseries[k][i])); fprintf(stderr, "\n");
		int j = 0;
		double temp = cabs(fseries[k][0]);
		for(i = 1; i < series[k]->data->length; i++) {
			if(temp < cabs(fseries[k][i])){
				j = i;
				temp = cabs(fseries[k][i]);
			}
		}
		fprintf(stderr, "signal %d %g\n", j, temp);

		j = 0;
		temp = cabs(fnseries[k][0]);
		for(i = 1; i < nseries[k]->length; i++) {
			if(temp < cabs(fnseries[k][i])){
				j = i;
				temp = cabs(fnseries[k][i]);
			}
		}
		fprintf(stderr, "noise  %d %g\n\n", j, temp);
#endif
	}


	/*
	 * Whitening
	 */


#if 0
//	double min = cabs(fnseries[0][0]);
//	for(k = 0; k < n; k++)
//		for(j = 0; j < (int) series[k]->data->length; j++)
//			if(min > cabs(fnseries[k][j]))
//				min = cabs(fnseries[k][j]);
//	min *= 1e3;
//	fprintf(stderr, "threshold = %g\n", min);

	for(k = 0; k < n; k++){
		double temp = 0;
		unsigned i;
		for(i = 0; i < series[k]->data->length; i++)
			temp += cabs(fnseries[k][i]);
		for(i = 0; i < series[k]->data->length; i++)
			fnseries[k][i] /= temp * temp;

		//threshold(fseries[k], fnseries[k], series[k]->data->length, min);
		whiten(fseries[k], fnseries[k], series[k]->data->length);
		//partial_whiten(fseries[k], fnseries[k], (int) (1024 * series[k]->data->length * series[k]->deltaT), (int) (2048 * series[k]->data->length * series[k]->deltaT));
		//filter(fseries[k], (int) (15 * series[k]->data->length * series[k]->deltaT), (int) (1024 * series[k]->data->length * series[k]->deltaT), series[k]->data->length);
	}
#endif
#if 0
	for(k = 0; k < n; k++){
		unsigned i; for(i = 0; i < series[k]->data->length; i++) fprintf(stderr, "%3d %g  %g\n", i, cabs(fseries[k][i]), cabs(fnseries[k][i])); fprintf(stderr, "\n");
		j = 0;
		double temp = cabs(fseries[k][0]);
		for(i = 1; i < series[k]->data->length; i++) {
			if(temp < cabs(fseries[k][i])){
				j = i;
				temp = cabs(fseries[k][i]);
			}
		}
		fprintf(stderr, "%d %g\n\n", j, temp);
	}
#endif


	/*
	 * Generate sky alm for each parameters of CBC
	 */


	*skyp = sh_series_new_zero(logprior->l_max, 0);
	if(generate_alm_sky(*skyp, fdplansp, fseries, fnseries, gmst_from_epoch_and_offset(series[0]->epoch, series[0]->data->length * series[0]->deltaT / 2.0), logprior)) {
		fprintf(stderr, "positive generate_alm_sky error\n");
		return -1;
	}


	*skyn = sh_series_new_zero(logprior->l_max, 0);
	if(generate_alm_sky(*skyn, fdplansn, fseries, fnseries, gmst_from_epoch_and_offset(series[0]->epoch, series[0]->data->length * series[0]->deltaT / 2.0), logprior)) {
		fprintf(stderr, "negative generate_alm_sky error\n");
		return -1;
	}


	/*
	 * Clean up
	 */


	for(k = 0; k < n; k++) {
		free(fseries[k]);
		free(fnseries[k]);
		fftw_destroy_plan(fftplans[k]);
		fftw_destroy_plan(nfftplans[k]);
	}

	return 0;
}
