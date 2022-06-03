/*
 * Copyright (C) 2022  Kipp C. Cannon
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
#include <math.h>
#include <stdio.h>

#include <sphradiometer/diagnostics.h>
#include <sphradiometer/sh_series.h>


static double min(double a, double b)
{
	return a < b ? a : b;
}


/*
 * ============================================================================
 *
 *                               Gradient Tests
 *
 * ============================================================================
 */


/*
 * the gradient of a band-limited impulse evaluated near the peak should
 * point towards the peak.
 */


static int grad_test_1(int l_max)
{
	int i;
	/* location of peak */
	double theta = acos(randrange(-1., +1.));
	double phi = randrange(0, 2. * M_PI);
	/* the impulse */
	struct sh_series *impulse = sh_series_impulse(l_max, theta, phi);
	/* the gradient */
	struct sh_series *u;
	struct sh_series *v;
	if(!sh_series_sintheta_grad(impulse, &u, &v)) {
		sh_series_free(impulse);
		return -1;
	}
	sh_series_free(impulse);

	for(i = 0; i < 1000; i++) {
		double x, y;
		/* a small displacement in some direction.  the physical
		 * displacement on sphere for displacements in phi is
		 * sin(theta) d phi, so d theta gets a factor of sin(theta)
		 * to keep the physical displacement in that direction in
		 * proportion to the d phi displacement */
		double alpha = randrange(-M_PI, M_PI);
		double d_theta = sin(theta) * 0.000001 * sin(alpha);
		double d_phi = 0.000001 * cos(alpha);
		/* don't let things get weird */
		if(theta + d_theta <= 0. || M_PI <= theta + d_theta ||
			phi + d_phi <= 0. || 2. * M_PI <= phi + d_phi)
			continue;

		/* evaluate gradient there */
		complex double u_val = sh_series_eval(u, theta + d_theta, phi + d_phi);
		complex double v_val = sh_series_eval(v, theta + d_theta, phi + d_phi);

		/* the impulse is real-valued, are the components of its
		 * gradient? */
		x = min(fabs(carg(u_val)), fabs(carg(u_val)) - M_PI);
		y = min(fabs(carg(v_val)), fabs(carg(v_val)) - M_PI);
		if(x > 1e-14 || y > 1e-14) {
			fprintf(stderr, "gradient not real valued: carg(u)=%.16g carg(v)=%.16g\n", x, y);
			goto fail;
		}
		/* clear the imaginary component for simplicity */
		u_val = creal(u_val);
		v_val = creal(v_val);

		/* because d_theta and sin(theta) d_phi are small, u and v
		 * should be proportional to them because for sufficiently
		 * small displacements spheres look flat.  therefore, we
		 * should be able to recover the direction in which we've
		 * been displaced from the peak from u and v. */

		x = atan2(-u_val, -v_val);
		if(fabs(alpha - x) > 1e-6) {
			fprintf(stderr, "gradient does not point towards peak:\n\tangle from +ve iso-latitude of displacement from peak %.16g ran\n\t-ve grad points in direction %.16g rad\n\t|error| = %.16g rad\n", alpha, x, fabs(alpha - x));
			goto fail;
		}
	}

	/* success */
	sh_series_free(u);
	sh_series_free(v);
	return 0;

fail:
	fprintf(stderr, "error occured testing gradient around impulse at theta = %.16g rad, phi = %.16g rad\n", theta, phi);
	sh_series_free(u);
	sh_series_free(v);
	return -1;
}


/*
 * compare the numerical derivatives to some results known analytically
 */


static double test_func(double theta, double phi, void *data)
{
	return pow(sin(theta), 2.) * sin(2. * phi);
}

static double test_func_u(double theta, double phi)
{
	/* sin theta d/d theta */
	return 2. * pow(sin(theta), 2.) * cos(theta) * sin(2. * phi);
}

static double test_func_v(double theta, double phi)
{
	/* d/d phi */
	return 2. * pow(sin(theta), 2.) * cos(2. * phi);
}


static int grad_test_2(void)
{
	struct sh_series *series = sh_series_from_realfunc(sh_series_new(13, 0), test_func, NULL);
	complex double *pixels;
	int ntheta, nphi;
	int i, j;
	double *cos_theta_array;
	double *cos_theta_weights;
	struct sh_series *u, *v;

	/* test mesh conversion just to be sure (confirms function has
	 * an exact representation in a finite Ylm basis) */
	pixels = sh_series_to_mesh(series);
	free(sh_series_mesh_new(series->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights));
	free(cos_theta_weights);
	for(j = 0; j < ntheta; j++) {
		double theta = acos(cos_theta_array[j]);
		for(i = 0; i < nphi; i++) {
			double phi = (2. * M_PI / nphi) * i;
			double correct = test_func(theta, phi, NULL);
			complex double p = pixels[j * nphi + i];
			double error = cabs(p - correct);
			if(error > 1e-14) {
				fprintf(stderr, "bad func at theta=%.16g phi=%.16g\n\tgot %.16g+I*%.16g expected %.16g |error|=%.16g\n", theta, phi, creal(p), cimag(p), correct, error);
				free(cos_theta_array);
				free(pixels);
				sh_series_free(series);
				return -1;
			}
		}
	}
	free(cos_theta_array);
	free(pixels);

	if(!sh_series_sintheta_grad(series, &u, &v)) {
		sh_series_free(series);
		return -1;
	}
	sh_series_free(series);

	/* derivative w.r.t. phi */
	pixels = sh_series_to_mesh(v);
	if(!pixels) {
		sh_series_free(u);
		sh_series_free(v);
		return -1;
	}
	free(sh_series_mesh_new(v->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights));
	free(cos_theta_weights);
	for(j = 0; j < ntheta; j++) {
		double theta = acos(cos_theta_array[j]);
		for(i = 0; i < nphi; i++) {
			double phi = (2. * M_PI / nphi) * i;
			double correct = test_func_v(theta, phi);
			complex double p = pixels[j * nphi + i];
			double error = cabs(p - correct) / fabs(correct);
			if(error > 1e-13) {
				fprintf(stderr, "bad d/dphi at theta=%.16g phi=%.16g\n\tgot %.16g+I*%.16g expected %.16g |error|=%.2g%%\n", theta, phi, creal(p), cimag(p), correct, 100. * error);
				free(cos_theta_array);
				free(pixels);
				sh_series_free(u);
				sh_series_free(v);
				return -1;
			}
		}
	}
	free(cos_theta_array);
	free(pixels);
	free(v);

	/* sin(theta) * deriviative w.r.t. theta */
	pixels = sh_series_to_mesh(u);
	free(sh_series_mesh_new(u->l_max, &ntheta, &nphi, &cos_theta_array, &cos_theta_weights));
	free(cos_theta_weights);
	for(j = 0; j < ntheta; j++) {
		double theta = acos(cos_theta_array[j]);
		for(i = 0; i < nphi; i++) {
			double phi = (2. * M_PI / nphi) * i;
			double correct = test_func_u(theta, phi);
			complex double p = pixels[j * nphi + i];
			double error = cabs(p - correct);
			if(error > 1e-14) {
				fprintf(stderr, "bad d/dtheta at theta=%.16g phi=%.16g\n\tgot %.16g+I*%.16g expected %.16g |error|=%.16g\n", theta, phi, creal(p), cimag(p), correct, error);
				free(cos_theta_array);
				free(pixels);
				sh_series_free(u);
				return -1;
			}
		}
	}
	free(cos_theta_array);
	free(pixels);
	free(u);

	return 0;
}


/*
 * ============================================================================
 *
 *                              Check Laplacian
 *
 * ============================================================================
 */


static double test_func_laplacian(double theta, double phi, void *data)
{
	/* 1/sin(theta) [ d/dtheta (sin(theta) d/dtheta) + 1/sin^2(theta) d^2/dphi^2 ] sin^2(theta) sin(2 phi)
	 *
	 * = 1/sin(theta) sin(2 phi) d/dtheta (sin(theta) 2 sin(theta)cos(theta)) - 4 sin(2 phi)
	 *
	 * = 1/sin(theta) 2 (-sin(theta) + 3 cos^2(theta) sin(theta)) sin(2 phi) - 4 sin(2 phi)
	 *
	 * = 2 [ (3 cos^2(theta) - 1) - 2 ] sin(2 phi)
	 *
	 * = 2 [ 3 cos^2(theta) - 3 ] sin(2 phi)
	 *
	 * = 6 [ cos^2(theta) - 1 ] sin(2 phi)
	 *
	 * = -6 sin^2(theta) sin(2 phi)
	 *
	 * = -6 * test_func() (it's a spherical harmonic, an eigenfunction
	 * of this operator)
	 */

	return -6. * test_func(theta, phi, data);
}


static int laplacian_test_1(void)
{
	struct sh_series *series = sh_series_from_realfunc(sh_series_new(13, 0), test_func, NULL);
	int i;

	if(!series)
		return -1;
	sh_series_laplacian(series);

	for(i = 0; i < 100; i++) {
		double theta = acos(randrange(-1., +1.));
		double phi = randrange(-M_PI, +M_PI);
		complex double correct = test_func_laplacian(theta, phi, NULL);
		complex double val = sh_series_eval(series, theta, phi);
		double error = cabs(correct - val);
		if(error > 1e-13) {
			fprintf(stderr, "laplacian failed at theta=%.16g phi=%.16g:\n\texpected %.16g+I*%.16g got %.16g+I*%.16g |error|=%.16g\n", theta, phi, creal(correct), cimag(correct), creal(val), cimag(val), error);
			sh_series_free(series);
			return -1;
		}
	}

	sh_series_free(series);
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
	int i;

	for(i = 0; i < 1000; i++)
		assert(grad_test_1(13) == 0);
	assert(grad_test_2() == 0);

	assert(laplacian_test_1() == 0);

	exit(0);
}
