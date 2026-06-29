/*
 * Copyright (C) 2026  Kipp C. Cannon
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
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <sphradiometer/sh_series.h>


/* FIXME: because of the repeated evaluations of Legendre polynomials, this
 * code would benefit from a multi-series evaluator API, something that
 * evaluates an array of functions at the same point returning an array of
 * their values.
 */


/*
 * ============================================================================
 *
 *                                  Helpers
 *
 * ============================================================================
 */


/*
 * the extremizer code doesn't know the rules for the co-ordinates and can
 * cross over wrap-around boundaries if the extremum is near one, so this
 * is used to sanitize co-ordinates before evaluating functions on the
 * sphere.  theta is put into [0, pi] and phi into [0, 2pi).  the values
 * are modified in-place.
 */


static void condition_theta_phi(double *theta, double *phi)
{
	/*fprintf(stderr, "condition_theta_phi(%.16g %.16g) --> ", *theta, *phi);*/
	/* unwind theta into [0, 2pi) */

	*theta = fmod(*theta, 2 * M_PI);
	if(*theta < 0.)
		*theta += 2 * M_PI;

	/* constrain theta to [0, pi], flipping phi to the other side of
	 * the sphere if needed */

	if(*theta > M_PI) {
		*theta = 2. * M_PI - *theta;
		*phi += M_PI;
	}

	/* finally unwind phi into [0, 2pi) */

	*phi = fmod(*phi, 2 * M_PI);
	if(*phi < 0.)
		*phi += 2 * M_PI;
	/*fprintf(stderr, "%.16g %.16g\n", *theta, *phi);*/

	assert(0. <= *theta && *theta <= M_PI);
	assert(0. <= *phi && *phi < 2. * M_PI);
}


/*
 * ============================================================================
 *
 *                   Wrapper Functions for Extremum Finding
 *
 * ============================================================================
 */


/*
 * params structure for use with GSL's "fdf" minimizers
 */


struct fdf_params {
	const struct sh_series *series;
	struct sh_series *u;
	struct sh_series *v;
};


/*
 * function and gradient evaluators to use GSL's "fdf" minimizers to find
 * the minimum
 */


static double f_for_min(const gsl_vector *x, void *params)
{
	struct fdf_params *fdf_params = params;
	double theta = gsl_vector_get(x, 0);
	double phi = gsl_vector_get(x, 1);

	condition_theta_phi(&theta, &phi);

	double val = creal(sh_series_eval(fdf_params->series, theta, phi));

	/*fprintf(stderr, "f_for_min(%.16g, %.16g) = %.16g\n", theta, phi, val);*/

	return val;
}


static void df_for_min(const gsl_vector *x, void *params, gsl_vector *g)
{
	struct fdf_params *fdf_params = params;
	double theta = gsl_vector_get(x, 0);
	double phi = gsl_vector_get(x, 1);

	condition_theta_phi(&theta, &phi);

	double sin_theta = sin(theta);

	double u = creal(sh_series_eval(fdf_params->u, theta, phi));
	double v = creal(sh_series_eval(fdf_params->v, theta, phi));

	/*fprintf(stderr, "grad f_for_min(%.16g, %.16g) = %.16g theta + %.16g phi\n", theta, phi, u / sin_theta, v / sin_theta);*/

	gsl_vector_set(g, 0, u / sin_theta);
	gsl_vector_set(g, 1, v / sin_theta);
}


static void fdf_for_min(const gsl_vector *x, void *params, double *f, gsl_vector *g)
{
	*f = f_for_min(x, params);
	df_for_min(x, params, g);
}


/*
 * function and gradient evaluators to use GSL's "fdf" minimizers to find
 * the maximum
 */


static double f_for_max(const gsl_vector *x, void *params)
{
	return -f_for_min(x, params);
}


static void df_for_max(const gsl_vector *x, void *params, gsl_vector *g)
{
	df_for_min(x, params, g);

	gsl_vector_set(g, 0, -gsl_vector_get(g, 0));
	gsl_vector_set(g, 1, -gsl_vector_get(g, 1));
}


static void fdf_for_max(const gsl_vector *x, void *params, double *f, gsl_vector *g)
{
	*f = f_for_max(x, params);
	df_for_max(x, params, g);
}


/*
 * find a point near the extremum
 */


static double find_near_max(const struct sh_series *series, double *theta, double *phi)
{
	struct sh_series *low_bandwidth;
	/* needed to work around const'edness.  optimizer should remove */
	const struct sh_series *const_low_bandwidth;
	int ntheta, nphi;
	int i, j;
	double *cos_theta_array;
	complex double *mesh, *m;
	double max;

	/*
	 * for safety
	 */

	*theta = *phi = NAN;

	/*
	 * make a low-bandwidth version of the function
	 */

	if(series->l_max > 6) {
		low_bandwidth = sh_series_copy(series);
		if(!low_bandwidth)
			return NAN;
		sh_series_resize(low_bandwidth, 6);
		if(!low_bandwidth)
			return NAN;
		const_low_bandwidth = low_bandwidth;
	} else
		const_low_bandwidth = series;

	/*
	 * evaluate the low-bandwidth version on a mesh
	 */

	mesh = sh_series_to_mesh(const_low_bandwidth, &ntheta, &nphi, &cos_theta_array);
	if(const_low_bandwidth != series)
		sh_series_free(low_bandwidth);
	if(!mesh)
		return NAN;

	/*
	 * select the maximum re(pixel)
	 */

	m = mesh;
	max = creal(*m);
	*theta = acos(cos_theta_array[0]);
	*phi = 0.;
	/*fprintf(stderr, "best guess max %.16g @ %.16g %.16g\n", max, *theta, *phi);*/
	for(j = 0; j < ntheta; j++) {
		for(i = 0; i < nphi; i++) {
			double val = creal(*m++);
			if(val > max) {
				*theta = acos(cos_theta_array[j]);
				*phi = i * 2 * M_PI / nphi;
				max = val;
				/*fprintf(stderr, "best guess max %.16g @ %.16g %.16g\n", max, *theta, *phi);*/
			}
		}
	}

	/*
	 * success
	 */

	free(mesh);
	free(cos_theta_array);

	return max;
}


/*
 * ============================================================================
 *
 *                                  Extrema
 *
 * ============================================================================
 */


/*
 * report the co-ordinates and value of the maximum of the real-valued
 * function on the sphere described by series.  on success the co-ordinates
 * are stored in the addresses pointed to by theta and phi, and the value
 * is returned.  on failure NaN is returned.
 *
 * NOTE:  this is a first attempt at this, it's not working very well, I'm
 * still experimenting.  use at your own risk.
 *
 * NOTE:  this uses GSL's implementation of the
 * Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm, which computes an
 * iterative approximation of the inverse of the function's Hessian matrix.
 * I believe we might be able to obtain a closed form expression for that
 * matrix, and if so then this code might be improved by writing a custom
 * gradient descent algorithm that takes advantage of the explicit form of
 * the inverse Hessian.
 */


double sh_series_real_maximum(const struct sh_series *series, double *theta, double *phi)
{
	struct fdf_params fdf_params = {
		.series = series,
	};
	gsl_multimin_function_fdf fdf = {
		.f = f_for_max,
		.df = df_for_max,
		.fdf = fdf_for_max,
		.n = 2,
		.params = &fdf_params
	};
	double val;
	gsl_vector *x = gsl_vector_alloc(2);
	if(!x)
		return NAN;

	/*
	 * find the co-ordinates of the maximum magnitude of a
	 * low-bandwidth approximation of series.  this will be the
	 * starting point for the gradient descent
	 */

	if(isnan(find_near_max(series, gsl_vector_ptr(x, 0), gsl_vector_ptr(x, 1)))) {
		gsl_vector_free(x);
		return NAN;
	}

	/*
	 * compute sin(theta) * grad of input series.  this completes the
	 * initialization of fdf_params.
	 */

	if(!sh_series_sintheta_grad(series, &fdf_params.u, &fdf_params.v)) {
		gsl_vector_free(x);
		return NAN;
	}

	/*
	 * initialize the "minimizer" (we're using it to find a maximum)
	 */

	gsl_multimin_fdfminimizer *extremizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, 2);
	if(!extremizer) {
		gsl_vector_free(x);
		sh_series_free(fdf_params.u);
		sh_series_free(fdf_params.v);
		return NAN;
	}

	gsl_multimin_fdfminimizer_set(extremizer, &fdf, x, M_PI / 20., 0.001);
	gsl_vector_free(x);

	for(int iter = 0; iter < 100; iter++) {
		int status;

		status = gsl_multimin_fdfminimizer_iterate(extremizer);
		if(status) {
			/*fprintf(stderr, "iterator failed\n");*/
			break;
		}

		status = gsl_multimin_test_gradient(extremizer->gradient, 1e-3);
		if(status == GSL_SUCCESS) {
			/*fprintf(stderr, "iterator converged\n");*/
			break;
		}

		/*x = gsl_multimin_fdfminimizer_x(extremizer);
		*theta = gsl_vector_get(x, 0);
		*phi = gsl_vector_get(x, 1);
		val = -gsl_multimin_fdfminimizer_minimum(extremizer);
		fprintf(stderr, "iter %d max %.16g @ %.16g %.16g\n", iter, val, *theta, *phi);*/
	}

	x = gsl_multimin_fdfminimizer_x(extremizer);
	*theta = gsl_vector_get(x, 0);
	*phi = gsl_vector_get(x, 1);
	val = -gsl_multimin_fdfminimizer_minimum(extremizer);

	gsl_multimin_fdfminimizer_free(extremizer);
	sh_series_free(fdf_params.u);
	sh_series_free(fdf_params.v);

	return val;
}
