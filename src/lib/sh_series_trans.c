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
#include <string.h>
#include <fftw3.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline2d.h>
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


static complex double cexpi(double x)
{
	/* about 1% faster than cexp(I * x) */
	return cos(x) + I * sin(x);
}


/*
 * Return Y_{lm}(theta, phi), and Y^{*}_{lm}(theta, phi).
 */


complex double sh_series_Y(unsigned int l, int m, double theta, double phi)
{
	assert(0. <= theta && theta <= M_PI);
	return ((m < 0 && m & 1) ? -1.0 : +1.0) * cexpi(m * phi) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
}


complex double sh_series_Yconj(unsigned int l, int m, double theta, double phi)
{
	assert(0. <= theta && theta <= M_PI);
	return ((m < 0 && m & 1) ? -1.0 : +1.0) * cexpi(-m * phi) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
}


/*
 * Compute an array of Y_{lm}(theta, phi), and Y^{*}_{lm}(theta, phi) for l
 * in [|m|, ..., l_max].
 */


complex double *sh_series_Y_array(complex double *array, unsigned int l_max, int m, double theta, double phi)
{
	int n = l_max - abs(m) + 1;
	double tmp[n];
	const complex double factor = ((m < 0) && (m & 1) ? -1.0 : +1.0) * cexpi(m * phi);
	int i;

	assert(0. <= theta && theta <= M_PI);
	gsl_sf_legendre_sphPlm_array(l_max, abs(m), cos(theta), tmp);

	for(i = 0; i < n; i++)
		array[i] = factor * tmp[i];

	return array;
}


complex double *sh_series_Yconj_array(complex double *array, unsigned int l_max, int m, double theta, double phi)
{
	int n = l_max + 1 - abs(m);
	double tmp[n];
	const complex double factor = ((m < 0) && (m & 1) ? -1.0 : +1.0) * cexpi(-m * phi);
	int i;

	assert(0. <= theta && theta <= M_PI);
	gsl_sf_legendre_sphPlm_array(l_max, abs(m), cos(theta), tmp);

	for(i = 0; i < n; i++)
		array[i] = factor * tmp[i];

	return array;
}


/*
 * ============================================================================
 *
 *                                 Pixel Mesh
 *
 * ============================================================================
 */


/*
 * Compute and return the parameters of the pixel mesh.  The pixel mesh is
 * rectangular in (theta, phi), with phi = i * (2 pi / nphi) for integer i
 * in [0, nphi), and cos(theta) being the ntheta = (2*l_max + 2)
 * Gauss-Legendre quadrature points on the interval [-1, +1].  The pixel
 * order is m[j * nphi + i] for co-ordinates (theta[j], phi[i]).
 */


static int pixels_from_l_max(unsigned l_max, int *ntheta, int *nphi, double **cos_theta_array, double **cos_theta_weights)
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
		return -1;
	}

	for(i = 0; i < *ntheta; i++)
		gsl_integration_glfixed_point(-1., +1., i, &(*cos_theta_array)[i], &(*cos_theta_weights)[i], t);

	gsl_integration_glfixed_table_free(t);

	/* reverse the cos_theta_array and cos_theta_weights so that the
	 * co-ordinates are ordered by theta instead of by cos(theta).  we
	 * don't really have to reverse the weights because they're
	 * symmetric about the origin */

	for(i = 0; i < *ntheta / 2; i++) {
		double tmp = (*cos_theta_array)[i];
		(*cos_theta_array)[i] = (*cos_theta_array)[*ntheta - 1 - i];
		(*cos_theta_array)[*ntheta - 1 - i] = tmp;
	}

	return 0;
}


static complex double *sh_series_mesh_new(unsigned int l_max, int *ntheta, int *nphi, double **cos_theta_array, double **cos_theta_weights)
{
	complex double *mesh;

	if(pixels_from_l_max(l_max, ntheta, nphi, cos_theta_array, cos_theta_weights))
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

	if(pixels_from_l_max(l_max, ntheta, nphi, cos_theta_array, cos_theta_weights))
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
 * ============================================================================
 *
 *                        sh_series Object Evaluation
 *
 * ============================================================================
 */


/*
 * Evaluate an sh_series object at theta, phi.
 */


complex double sh_series_eval(const struct sh_series *series, double theta, double phi)
{
	complex double vals[series->l_max + 1];
	complex double *v_last = vals + series->l_max + 1;
	complex double *coeff = series->coeff;
	int m_max = series->polar ? 0 : series->l_max;
	int m;
	complex double val = 0.0;

	for(m = 0; m <= m_max; m++, v_last--) {
		complex double *v;
		complex double x = 0.;
		sh_series_Y_array(vals, series->l_max, m, theta, phi);
		for(v = vals; v < v_last;)
			x += *(coeff++) * *(v++);
		val += x;
	}
	for(m = -m_max, v_last++; m < 0; m++, v_last++) {
		complex double *v;
		complex double x = 0.;
		sh_series_Y_array(vals, series->l_max, m, theta, phi);
		for(v = vals; v < v_last;)
			x += *(coeff++) * *(v++);
		val += x;
	}

	return val;
}


/*
 * sh_series_eval_interp
 */


struct sh_series_eval_interp {
	int ntheta, nphi;
	double *theta;
	gsl_interp_accel *theta_acc;
	double *phi;
	gsl_interp_accel *phi_acc;
	gsl_spline2d *re;
	gsl_spline2d *im;
};


struct sh_series_eval_interp *sh_series_eval_interp_new(const struct sh_series *series)
{
	struct sh_series_eval_interp *interp = malloc(sizeof(*interp));
	double *tmp;
	complex double *mesh;
	double *interp_mesh;
	int i, j;
	if(!interp)
		return NULL;

	if(pixels_from_l_max(series->l_max, &interp->ntheta, &interp->nphi, &interp->theta, &tmp)) {
		free(interp);
		return NULL;
	}
	free(tmp);

	for(j = 0; j < interp->ntheta; j++)
		interp->theta[j] = acos(interp->theta[j]);

	/* to ensure periodicity in phi, 4 extra columns are added and two
	 * rows.  the real array of phi is the central nphi columns, with
	 * two columns on the left and on the right providing duplicate
	 * copies of the other end of the array thereby inducing the
	 * interpolator to construct a periodic function.  the theta array
	 * gets the north and south poles added to prevent "y value out of
	 * range" errors in the GSL interpolator. */
	tmp = realloc(interp->theta, (interp->ntheta + 2) * sizeof(*interp->theta));
	interp->phi = malloc((interp->nphi + 4) * sizeof(*interp->phi));
	interp->re = gsl_spline2d_alloc(gsl_interp2d_bicubic, interp->nphi + 4, interp->ntheta + 2);
	interp->im = gsl_spline2d_alloc(gsl_interp2d_bicubic, interp->nphi + 4, interp->ntheta + 2);
	interp->theta_acc = gsl_interp_accel_alloc();
	interp->phi_acc = gsl_interp_accel_alloc();
	mesh = sh_series_to_mesh(series);
	interp_mesh = malloc((interp->nphi + 4) * (interp->ntheta + 2) * sizeof(*interp_mesh));
	if(!tmp || !interp->phi || !interp->re || !interp->im || !interp->theta_acc || !interp->phi_acc || !mesh || !interp_mesh) {
		free(interp->theta);
		free(interp->phi);
		gsl_spline2d_free(interp->re);
		gsl_spline2d_free(interp->im);
		gsl_interp_accel_free(interp->theta_acc);
		gsl_interp_accel_free(interp->phi_acc);
		free(mesh);
		free(interp_mesh);
		free(interp);
		return NULL;
	}
	interp->theta = tmp;

	/* === can now use sh_series_eval_interp_free() for clean up === */

	/* add poles to theta co-ordinate array */

	memmove(interp->theta + 1, interp->theta, interp->ntheta * sizeof(*interp->theta));
	interp->theta[0] = 0.;
	interp->theta[interp->ntheta + 1] = M_PI;

	/* initialize the phi co-ordinate array */
	for(i = -2; i < interp->nphi + 2; i++)
		interp->phi[i + 2] = (2. * M_PI / interp->nphi) * i;

	/* initialize the real part interpolator */
	{
	double x;
	x = creal(sh_series_eval(series, 0., 0.));
	for(i = 0; i < interp->nphi + 4; i++)
		interp_mesh[(0) * (interp->nphi + 4) + i] = x;
	x = creal(sh_series_eval(series, M_PI, 0.));
	for(i = 0; i < interp->nphi + 4; i++)
		interp_mesh[(interp->ntheta + 1) * (interp->nphi + 4) + i] = x;
	}
	for(j = 0; j < interp->ntheta; j++) {
		interp_mesh[(j + 1) * (interp->nphi + 4) + 0] = creal(mesh[j * interp->nphi + interp->nphi - 2]);
		interp_mesh[(j + 1) * (interp->nphi + 4) + 1] = creal(mesh[j * interp->nphi + interp->nphi - 1]);
		for(i = 0; i < interp->nphi; i++)
			interp_mesh[(j + 1) * (interp->nphi + 4) + i + 2] = creal(mesh[j * interp->nphi + i]);
		interp_mesh[(j + 1) * (interp->nphi + 4) + interp->nphi + 2] = creal(mesh[j * interp->nphi + 0]);
		interp_mesh[(j + 1) * (interp->nphi + 4) + interp->nphi + 3] = creal(mesh[j * interp->nphi + 1]);
		}
	if(gsl_spline2d_init(interp->re, interp->phi, interp->theta, interp_mesh, interp->nphi + 4, interp->ntheta + 2)) {
		free(mesh);
		free(interp_mesh);
		sh_series_eval_interp_free(interp);
		return NULL;
	}

	/* initialize the imaginary part interpolator */
	{
	double x;
	x = creal(sh_series_eval(series, 0., 0.));
	for(i = 0; i < interp->nphi + 4; i++)
		interp_mesh[(0) * (interp->nphi + 4) + i] = x;
	x = cimag(sh_series_eval(series, M_PI, 0.));
	for(i = 0; i < interp->nphi + 4; i++)
		interp_mesh[(interp->ntheta + 1) * (interp->nphi + 4) + i] = x;
	}
	for(j = 0; j < interp->ntheta; j++) {
		interp_mesh[(j + 1) * (interp->nphi + 4) + 0] = cimag(mesh[j * interp->nphi + interp->nphi - 2]);
		interp_mesh[(j + 1) * (interp->nphi + 4) + 1] = cimag(mesh[j * interp->nphi + interp->nphi - 1]);
		for(i = 0; i < interp->nphi; i++)
			interp_mesh[(j + 1) * (interp->nphi + 4) + i + 2] = cimag(mesh[j * interp->nphi + i]);
		interp_mesh[(j + 1) * (interp->nphi + 4) + interp->nphi + 2] = cimag(mesh[j * interp->nphi + 0]);
		interp_mesh[(j + 1) * (interp->nphi + 4) + interp->nphi + 3] = cimag(mesh[j * interp->nphi + 1]);
		}
	if(gsl_spline2d_init(interp->im, interp->phi, interp->theta, interp_mesh, interp->nphi + 4, interp->ntheta + 2)) {
		free(mesh);
		free(interp_mesh);
		sh_series_eval_interp_free(interp);
		return NULL;
	}

	free(mesh);
	free(interp_mesh);

	return interp;
}


void sh_series_eval_interp_free(struct sh_series_eval_interp *interp)
{
	if(interp) {
		free(interp->theta);
		free(interp->phi);
		gsl_spline2d_free(interp->re);
		gsl_spline2d_free(interp->im);
		gsl_interp_accel_free(interp->theta_acc);
		gsl_interp_accel_free(interp->phi_acc);
	}
	free(interp);
}


double complex sh_series_eval_interp(const struct sh_series_eval_interp *interp, double theta, double phi)
{
	if(phi < 0.)
		phi += 2. * M_PI * (floor(phi / (2. * M_PI)) + 1.);
	else if(phi >= 2. * M_PI)
		phi -= 2. * M_PI * floor(phi / (2. * M_PI));
	assert(0. <= theta && theta <= M_PI);
	return gsl_spline2d_eval(interp->re, phi, theta, interp->phi_acc, interp->theta_acc) + I * gsl_spline2d_eval(interp->im, phi, theta, interp->phi_acc, interp->theta_acc);
}


/*
 * ============================================================================
 *
 *                       Spherical Harmonic Transforms
 *
 * ============================================================================
 */


static complex double *sh_series_mesh_from_func(unsigned int l_max, complex double (*func)(double, double, void *), void *data, int *ntheta, int *nphi)
{
	int nt, np;
	double *cos_theta_array, *cos_theta_weights;
	complex double *f = sh_series_mesh_new(l_max, &nt, &np, &cos_theta_array, &cos_theta_weights);
	const double dphi = 2. * M_PI / np;

	if(!f)
		return NULL;

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	{
	complex double *_f = f;
	for(int i = 0; i < nt; i++) {
		double theta = acos(cos_theta_array[i]);
		for(int j = 0; j < np; j++)
			*(_f++) = func(theta, dphi * j, data);
	}
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
	const double dphi = 2. * M_PI / np;

	if(!f)
		return NULL;

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	{
	double *_f = f;
	for(int i = 0; i < nt; i++) {
		double theta = acos(cos_theta_array[i]);
		for(int j = 0; j < np; j++)
			*(_f++) = func(theta, dphi * j, data);
	}
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
	unsigned int bw = series->l_max + 1;
	int m_max = series->polar ? 0 : series->l_max;
	complex double *F;
	double *P;

	if(pixels_from_l_max(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights))
		return NULL;
	F = malloc(ntheta * nphi * sizeof(*F));
	P = malloc(ntheta * (series->l_max + 1) * sizeof(*P));
	if(!F || !P) {
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

	{
	complex double *_F = F;
	for(int i = 0; i < ntheta; i++) {
		double weight = cos_theta_weights[i] * (2. * M_PI / nphi);
		for(int j = 0; j < nphi; j++)
			*(_F++) *= weight;
	}
	}

	/* for each non-negative m, */
	for(int m = 0; m <= m_max; m++) {
		/* populate P[][] with (normalized) P_{l m}(cos theta) */
		for(int i = 0; i < ntheta; i++)
			gsl_sf_legendre_sphPlm_array(series->l_max, m, cos_theta_array[i], P + i * bw + m);
		/* for each l s.t. m <= l <= l_max, */
		for(int l = m; l <= (int) series->l_max; l++) {
			/* integrate H_{m}(theta) sin(theta) P_{l m}(cos
			 * theta) over theta (+m case) */
			complex double *_F = F + m;
			double *_P = P + l;
			complex double c = 0.0;
			for(int i = 0; i < ntheta; i++, _F += nphi, _P += bw)
				c += *_F * *_P;
			sh_series_set(series, l, m, c);
		}
		if(m) {
			const double sign = m & 1 ? -1.0 : +1.0;
			/* for each l s.t. m <= l <= l_max, */
			for(int l = m; l <= (int) series->l_max; l++) {
				/* integrate H_{-m}(theta) sin(theta) P_{l
				 * -m}(cos theta) over theta (-m case) */
				complex double *_F = F + nphi - m;
				double *_P = P + l;
				complex double c = 0.0;
				for(int i = 0; i < ntheta; i++, _F += nphi, _P += bw)
					c += *_F * *_P;
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
	unsigned int bw = series->l_max + 1;
	int m_max = series->polar ? 0 : series->l_max;
	complex double *F;
	double *P;

	if(pixels_from_l_max(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights))
		return NULL;
	F = malloc(ntheta * (nphi / 2 + 1) * sizeof(*F));
	P = malloc(ntheta * (series->l_max + 1) * sizeof(*P));
	if(!F || !P) {
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

	{
	complex double *_F = F;
	for(int i = 0; i < ntheta; i++) {
		double weight = cos_theta_weights[i] * (2. * M_PI / nphi);
		for(int j = 0; j < (nphi / 2 + 1); j++)
			*(_F++) *= weight;
	}
	}

	/* for each non-negative m, */
	for(int m = 0; m <= m_max; m++) {
		/* populate P[][] with (normalized) P_{lm}(cos theta) */
		for(int i = 0; i < ntheta; i++)
			gsl_sf_legendre_sphPlm_array(series->l_max, m, cos_theta_array[i], P + i * bw + m);
		/* for each l s.t. m <= l <= l_max, */
		for(int l = m; l <= (int) series->l_max; l++) {
			/* integrate H_{m}(theta) sin(theta) P_{lm}(cos
			 * theta) over theta for +/-m */
			complex double *_F = F + m;
			double *_P = P + l;
			complex double c = 0.0;
			for(int i = 0; i < ntheta; i++, _F += nphi / 2 + 1, _P += bw)
				c += *_F * *_P;
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
	int m_max = series->polar ? 0 : series->l_max;
	complex double *mesh = sh_series_mesh_new(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights);
	complex double *F = calloc(ntheta * nphi, sizeof(*F));
	double P[series->l_max + 1];

	/* not needed */
	free(cos_theta_weights);

	if(!mesh || !F) {
		free(mesh);
		free(F);
		free(cos_theta_array);
		return NULL;
	}

	/* populate F[] */
	for(int i = 0; i < ntheta; i++)
		for(int m = -m_max; m <= m_max; m++) {
			complex double *coeff = series->coeff + sh_series_moffset(series->l_max, m);
			complex double x = 0.;
			/* compute (normalized) P_{lm}(cos theta) */
			gsl_sf_legendre_sphPlm_array(series->l_max, abs(m), cos_theta_array[i], P + abs(m));
			for(int l = abs(m); l <= (int) series->l_max; l++)
				x += *coeff++ * P[l];
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
 * l_max.  the peak of the impulse is located at the co-ordinates (theta,
 * phi).  returns a newly allocated sh_series on success, or NULL on
 * failure.  this is the spherical equivalent of the sinc function for 1-D
 * time series analysis, and is obtained in the same manner by projecting a
 * Dirac delta function onto the harmonic basis and retaining only those
 * components at or below the requested l_max.
 */


struct sh_series *sh_series_impulse(unsigned int l_max, double theta, double phi)
{
	struct sh_series *series = sh_series_new(l_max, 0);

	if(!series)
		return NULL;

	/* project a Dirac delta onto the spherical harmonic basis upto and
	 * including order l_max */

	for(int m = -(int) l_max; m <= (int) l_max; m++)
		sh_series_Yconj_array(series->coeff + sh_series_moffset(l_max, m), l_max, m, theta, phi);

	return series;
}
