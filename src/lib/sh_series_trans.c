/*
 * Copyright (C) 2006--2009,2019  Kipp C. Cannon
 * Copyright (C) 2019  Takuya Tsutsui
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
#include <stdlib.h>
#include <fftw3.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                       Spherical Harmonic Evaluation
 *
 * ============================================================================
 */


/*
 * Return Y_{lm}(theta, phi), and Y^{*}_{lm}(theta, phi).
 */


complex double sh_series_Y(unsigned int l, int m, double theta, double phi)
{
	assert(0. <= theta && theta <= M_PI);
	return ((m < 0 && m & 1) ? -1.0 : +1.0) * (cos(m * phi) + I * sin(m * phi)) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
}


complex double sh_series_Yconj(unsigned int l, int m, double theta, double phi)
{
	assert(0. <= theta && theta <= M_PI);
	return ((m < 0 && m & 1) ? -1.0 : +1.0) * (cos(m * phi) - I * sin(m * phi)) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
}


/*
 * Compute an array of Y_{lm}(theta, phi), and Y^{*}_{lm}(theta, phi) for l
 * in [|m|, ..., l_max].
 */


complex double *sh_series_Y_array(complex double *array, unsigned int l_max, int m, double theta, double phi)
{
	double tmp[l_max + 1 - abs(m)];
	const complex double factor = ((m < 0) && (m & 1) ? -1.0 : +1.0) * cexp(I * m * phi);
	int i;

	assert(0. <= theta && theta <= M_PI);
	gsl_sf_legendre_sphPlm_array(l_max, abs(m), cos(theta), tmp);

	for(i = 0; i <= (int) l_max - abs(m); i++)
		array[i] = factor * tmp[i];

	return array;
}


complex double *sh_series_Yconj_array(complex double *array, unsigned int l_max, int m, double theta, double phi)
{
	double tmp[l_max + 1 - abs(m)];
	const complex double factor = ((m < 0) && (m & 1) ? -1.0 : +1.0) * cexp(I * -m * phi);
	int i;

	assert(0. <= theta && theta <= M_PI);
	gsl_sf_legendre_sphPlm_array(l_max, abs(m), cos(theta), tmp);

	for(i = 0; i <= (int) l_max - abs(m); i++)
		array[i] = factor * tmp[i];

	return array;
}


/*
 * ============================================================================
 *
 *                       Spherical Harmonic Transforms
 *
 * ============================================================================
 */


/*
 * Evaluate an sh_series object at theta, phi.
 */


complex double sh_series_eval(const struct sh_series *series, double theta, double phi)
{
	complex double vals[series->l_max + 1];
	complex double *coeff = series->coeff;
	int l, m;
	complex double val = 0.0;

	sh_series_Y_array(vals, series->l_max, 0, theta, phi);
	for(l = 0; l <= (int) series->l_max; l++)
		val += *(coeff++) * vals[l];
	if(!series->polar) {
		for(m = 1; m <= (int) series->l_max; m++) {
			complex double *v = vals;
			sh_series_Y_array(v, series->l_max, m, theta, phi);
			for(l = m; l <= (int) series->l_max; l++)
				val += *(coeff++) * *(v++);
		}
		for(m = -(int) series->l_max; m < 0; m++) {
			complex double *v = vals;
			sh_series_Y_array(v, series->l_max, m, theta, phi);
			for(l = -m; l <= (int) series->l_max; l++)
				val += *(coeff++) * *(v++);
		}
	}

	return val;
}


/*
 * Given a function of theta, phi, evaluate it on a 2-dimensional mesh of
 * points, for later input to a projection function.
 */


static void samples_from_l_max(unsigned l_max, int *ntheta, int *nphi, double **cos_theta_array, double **cos_theta_weights)
{
	unsigned int bw = l_max + 1;

	gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(2 * bw);
	int i;

	*ntheta = 2 * bw;
	*nphi = 2 * bw;

	*cos_theta_array = malloc(*ntheta * sizeof(**cos_theta_array));
	*cos_theta_weights = malloc(*ntheta * sizeof(**cos_theta_weights));
	if(!t || !*cos_theta_array || !*cos_theta_weights) {
		gsl_integration_glfixed_table_free(t);
		free(*cos_theta_array);
		free(*cos_theta_weights);
		*cos_theta_array = *cos_theta_weights = NULL;
		return;
	}

	for(i = 0; i < *ntheta; i++)
		gsl_integration_glfixed_point(-1., +1., i, &(*cos_theta_array)[i], &(*cos_theta_weights)[i], t);

	gsl_integration_glfixed_table_free(t);
}


static complex double *sh_series_mesh_new(unsigned int l_max, int *ntheta, int *nphi, double **cos_theta_array, double **cos_theta_weights)
{
	complex double *mesh;

	samples_from_l_max(l_max, ntheta, nphi, cos_theta_array, cos_theta_weights);
	if(!*cos_theta_array)
		return NULL;

	mesh = malloc(*ntheta * *nphi * sizeof(*mesh));
	if(mesh)
		return mesh;

	free(*cos_theta_array);
	free(*cos_theta_weights);
	free(mesh);
	*cos_theta_array = *cos_theta_weights = NULL;
	return NULL;
}


static double *sh_series_real_mesh_new(unsigned int l_max, int *ntheta, int *nphi, double **cos_theta_array, double **cos_theta_weights)
{
	double *mesh;

	samples_from_l_max(l_max, ntheta, nphi, cos_theta_array, cos_theta_weights);
	if(!*cos_theta_array)
		return NULL;

	mesh = malloc(*ntheta * *nphi * sizeof(*mesh));
	if(mesh)
		return mesh;

	free(*cos_theta_array);
	free(*cos_theta_weights);
	free(mesh);
	*cos_theta_array = *cos_theta_weights = NULL;
	return NULL;
}


/*
 * This fills the array with samples of f in the same order and evaluated
 * at the same points on the sphere as required by the functions in S2kit.
 */


static complex double *sh_series_mesh_from_func(unsigned int l_max, complex double (*func)(double, double, void *), void *data, int *ntheta, int *nphi)
{
	int nt, np;
	double *cos_theta_array, *cos_theta_weights;
	complex double *f = sh_series_mesh_new(l_max, &nt, &np, &cos_theta_array, &cos_theta_weights);
	const double dphi = 2 * M_PI / np;
	int i, j;

	if(!f)
		return NULL;

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	for(i = 0; i < nt; i++) {
		double theta = acos(cos_theta_array[i]);
		for(j = 0; j < np; j++)
			*(f + i * np + j) = func(theta, dphi * j, data);
	}

	free(cos_theta_array);
	free(cos_theta_weights);

	return f;
}


static double *sh_series_mesh_from_realfunc(unsigned int l_max, double (*func)(double, double, void *), void *data, int *ntheta, int *nphi)
{
	int nt, np;
	double *cos_theta_array, *cos_theta_weights;
	double *f = sh_series_real_mesh_new(l_max, &nt, &np, &cos_theta_array, &cos_theta_weights);
	const double dphi = 2 * M_PI / np;
	int i, j;

	if(!f)
		return NULL;

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	for(i = 0; i < nt; i++) {
		double theta = acos(cos_theta_array[i]);
		for(j = 0; j < np; j++)
			*(f + i * np + j) = func(theta, dphi * j, data);
	}

	free(cos_theta_array);
	free(cos_theta_weights);

	return f;
}


/*
 * Given a 2-D mesh of samples, project the data onto spherical harmonics,
 * and return the coefficients as an sh_series object.  Returns NULL on
 * failure.
 */


struct sh_series *sh_series_from_mesh(struct sh_series *series, complex double *mesh)
{
	int ntheta, nphi;
	double *cos_theta_array, *cos_theta_weights;
	complex double *F;
	double *P;
	int l, m, i, j;

	samples_from_l_max(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights);
	F = malloc(ntheta * nphi * sizeof(*F));
	P = malloc(ntheta * (series->l_max + 1) * sizeof(*P));
	if(!F || !P || !cos_theta_array) {
		free(F);
		free(P);
		free(cos_theta_array);
		free(cos_theta_weights);
		return NULL;
	}

	/* Fourier transform the phi co-ordinate.  after this F[][]
	 * contains H_{m}(theta) */

	{
	int n[] = {nphi};
	fftw_plan plan = fftw_plan_many_dft(1, n, ntheta, mesh, NULL, 1, nphi, F, NULL, 1, nphi, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	}

	/* multiply the samples by dphi to complete the normalization of
	 * the Fourier transform, and also by the Gauss-Legendre weights
	 * which are needed for the integral over cos(theta) below so we
	 * don't do that in the inner loop */

	for(i = 0; i < ntheta; i++)
		for(j = 0; j < nphi; j++)
			*(F + i * nphi + j) *= cos_theta_weights[i] * (2 * M_PI / nphi);

	/* for each non-negative m, */
	for(m = 0; m <= (int) (series->polar ? 0 : series->l_max); m++) {
		/* populate P[][] with (normalized) P_{l m}(cos theta) */
		for(i = 0; i < ntheta; i++)
			gsl_sf_legendre_sphPlm_array(series->l_max, m, cos_theta_array[i], P + i * (series->l_max + 1) + m);
		/* for each l s.t. m <= l <= l_max, */
		for(l = m; l <= (int) series->l_max; l++) {
			/* integrate H_{m}(theta) sin(theta) P_{l m}(cos
			 * theta) over theta (+m case) */
			complex double c = 0.0;
			for(i = 0; i < ntheta; i++)
				c += *(F + i * nphi + m) * *(P + i * (series->l_max + 1) + l);
			sh_series_set(series, l, m, c);
		}
		if(m) {
			const double sign = m & 1 ? -1.0 : +1.0;
			/* for each l s.t. m <= l <= l_max, */
			for(l = m; l <= (int) series->l_max; l++) {
				/* integrate H_{-m}(theta) sin(theta) P_{l
				 * -m}(cos theta) over theta (-m case) */
				complex double c = 0.0;
				for(i = 0; i < ntheta; i++)
					c += *(F + i * nphi + (nphi - m)) * *(P + i * (series->l_max + 1) + l);
				sh_series_set(series, l, -m, sign * c);
			}
		}
	}

	/* clean up */
	free(P);
	free(F);
	free(cos_theta_array);
	free(cos_theta_weights);

	return series;
}


struct sh_series *sh_series_from_realmesh(struct sh_series *series, double *mesh)
{
	int ntheta, nphi;
	double *cos_theta_array, *cos_theta_weights;
	complex double *F;
	double *P;
	int l, m, i, j;

	samples_from_l_max(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights);
	F = malloc(ntheta * (nphi / 2 + 1) * sizeof(*F));
	P = malloc(ntheta * (series->l_max + 1) * sizeof(*P));
	if(!F || !P || !cos_theta_array) {
		free(F);
		free(P);
		free(cos_theta_array);
		free(cos_theta_weights);
		return NULL;
	}

	/* Fourier transform the phi co-ordinate.  after this F[][]
	 * contains H_{m}(theta) */

	{
	int n[] = {nphi};
	fftw_plan plan = fftw_plan_many_dft_r2c(1, n, ntheta, mesh, NULL, 1, nphi, F, NULL, 1, nphi / 2 + 1, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	}

	/* multiply the samples by dphi to complete the normalization of
	 * the Fourier transform, and also by the Gauss-Legendre weights
	 * which are needed for the integral over cos(theta) below so we
	 * don't do that in the inner loop */

	for(i = 0; i < ntheta; i++)
		for(j = 0; j < (nphi / 2 + 1); j++)
			*(F + i * (nphi / 2 + 1) + j) *= cos_theta_weights[i] * (2 * M_PI / nphi);

	/* for each non-negative m, */
	for(m = 0; m <= (int) (series->polar ? 0 : series->l_max); m++) {
		/* populate P[][] with (normalized) P_{lm}(cos theta) */
		for(i = 0; i < ntheta; i++)
			gsl_sf_legendre_sphPlm_array(series->l_max, m, cos_theta_array[i], P + i * (series->l_max + 1) + m);
		/* for each l s.t. m <= l <= l_max, */
		for(l = m; l <= (int) series->l_max; l++) {
			/* integrate H_{m}(theta) sin(theta) P_{lm}(cos
			 * theta) over theta for +/-m */
			complex double c = 0.0;
			for(i = 0; i < ntheta; i++)
				c += *(F + i * (nphi / 2 + 1) + m) * *(P + i * (series->l_max + 1) + l);
			sh_series_set(series, l, -m, (m & 1 ? -1 : 1) * conj(sh_series_set(series, l, m, c)));
		}
	}

	/* clean up */
	free(P);
	free(F);
	free(cos_theta_array);
	free(cos_theta_weights);

	return series;
}


/*
 * Given a function of theta, phi, project the function onto spherical
 * harmonics upto and including l = l_max, and return the coefficients as
 * an sh_series object.  Returns NULL on failure.
 */


struct sh_series *sh_series_from_func(struct sh_series *series, complex double (*func)(double, double, void *), void *data)
{
	complex double *mesh = sh_series_mesh_from_func(series->l_max, func, data, NULL, NULL);

	if(!mesh)
		return NULL;

	series = sh_series_from_mesh(series, mesh);
	free(mesh);

	return series;
}


struct sh_series *sh_series_from_realfunc(struct sh_series *series, double (*func)(double, double, void *), void *data)
{
	double *mesh = sh_series_mesh_from_realfunc(series->l_max, func, data, NULL, NULL);

	if(!mesh)
		return NULL;

	series = sh_series_from_realmesh(series, mesh);
	free(mesh);
	return series;
}


/*
 * Compute a pixel mesh from an sh_series object.  The inverse of
 * sh_series_from_mesh().
 */


complex double *sh_series_to_mesh(const struct sh_series *series)
{
	int ntheta, nphi;
	double *cos_theta_array, *cos_theta_weights;
	complex double *mesh = sh_series_mesh_new(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights);
	complex double *F = calloc(ntheta * nphi, sizeof(*F));
	double P[series->l_max + 1];
	int i;
	int l, m;

	/* not needed */
	free(cos_theta_weights);

	if(!mesh || !F) {
		free(mesh);
		free(F);
		free(cos_theta_array);
		return NULL;
	}

	/* populate F[] */
	for(i = 0; i < ntheta; i++)
		for(m = (series->polar ? 0 : -(int) series->l_max); m <= (series->polar ? 0 : (int) series->l_max); m++) {
			complex double *series_vals = series->coeff + sh_series_moffset(series->l_max, m);
			complex double x;
			/* compute (normalized) P_{lm}(cos theta) */
			gsl_sf_legendre_sphPlm_array(series->l_max, abs(m), cos_theta_array[i], P + abs(m));
			x = 0.;
			for(l = abs(m); l <= (int) series->l_max; l++)
				x += *series_vals++ * P[l];
			*(F + i * nphi + (m + nphi) % nphi) = (m < 0) && (m & 1) ? -x : +x;
		}

	/* Fourier transform the m co-ordinate to the phi co-ordinate */
	{
	int n[] = {nphi};
	fftw_plan plan = fftw_plan_many_dft(1, n, ntheta, F, NULL, 1, nphi, mesh, NULL, 1, nphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	}

	/* clean up */
	free(F);
	free(cos_theta_array);

	return mesh;
}


/*
 * generate a band-limited impulse on the sphere, with bandwidth set by
 * l_max.  the impulse is located at the co-ordinates (theta, phi).
 * returns a newly allocated sh_series on success, or NULL on failure.
 */


struct sh_series *sh_series_impulse(unsigned int l_max, double theta, double phi)
{
	struct sh_series *series = sh_series_new(l_max, 0);
	int m;

	if(!series)
		return NULL;

	/* project a Dirac delta onto the spherical harmonic basis upto and
	 * including order l_max */

	for(m = -(int) l_max; m <= (int) l_max; m++)
		sh_series_Yconj_array(series->coeff + sh_series_moffset(l_max, m), l_max, m, theta, phi);

	return series;
}
