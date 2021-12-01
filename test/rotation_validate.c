/*
 * Copyright (C) 2019  Kipp C. Cannon
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
#include <stdlib.h>
#include <string.h>

#include <sphradiometer/diagnostics.h>
#include <sphradiometer/sky.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                        Euler Rotation Matrix Tests
 *
 * ============================================================================
 */


#if 0
static double abs_fractional_difference(double x, double y)
{
	x = fabs(x);
	y = fabs(y);
	return (x == 0. && y == 0.) ?  0. : fabs(x - y) / (x > y ? x : y);
}
#endif


static void R_print(const double *R)
{
	fprintf(stderr, "%.5g %.5g %.5g\n", R[0], R[1], R[2]);
	fprintf(stderr, "%.5g %.5g %.5g\n", R[3], R[4], R[5]);
	fprintf(stderr, "%.5g %.5g %.5g\n", R[6], R[7], R[8]);
}


static int R_cmp(const double *R1, const double *R2, double tolerance)
{
	double max_delta = fabs(R1[0] - R2[0]);
	int i;
	for(i = 1; i < 9; i++) {
		double delta = fabs(R1[i] - R2[i]);
		if(delta > max_delta)
			max_delta = delta;
	}

	if(max_delta > tolerance) {
		fprintf(stderr, "R1 =\n");
		R_print(R1);
		fprintf(stderr, "!= R2 =\n");
		R_print(R2);
		fprintf(stderr, "max delta = %g, tolerance was %g\n", max_delta, tolerance);
	}

	return max_delta > tolerance;
}


static double *R_mult(double *dst, const double *src)
{
	double tmp[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	int i, j, k;
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			for(k = 0; k < 3; k++)
				tmp[3 * i + j] += dst[3 * i + k] * src[3 * k + j];
	memcpy(dst, tmp, sizeof(tmp));
	return dst;
}


static double *R_ident(void)
{
	/* identity */
	static double R[] = {
		1., 0., 0.,
		0., 1., 0.,
		0., 0., 1.
	};
	return R;
}


static double *R_gamma_eq_pi2(void)
{
	/* rotation by pi/2 about z axis:  x-->y, y-->-x */
	static double R[] = {
		0., -1., 0.,
		1.,  0., 0.,
		0.,  0., 1.
	};
	return R;
}


static double *R_beta_eq_pi2(void)
{
	/* rotation by pi/2 about y axis:  z-->x, x-->-z */
	static double R[] = {
		 0., 0., 1.,
		 0., 1., 0.,
		-1., 0., 0.
	};
	return R;
}


static double *R_alpha_eq_pi2(void)
{
	return R_gamma_eq_pi2();
}


static int test_euler_rotation_matrix(void)
{
	double *R;

	R = euler_rotation_matrix(0., 0., M_PI_2);
	assert(R_cmp(R_alpha_eq_pi2(), R, 1e-16) == 0);
	free(R);

	R = euler_rotation_matrix(0., M_PI_2, 0.);
	assert(R_cmp(R_beta_eq_pi2(), R, 1e-16) == 0);
	free(R);

	R = euler_rotation_matrix(M_PI_2, 0., 0.);
	assert(R_cmp(R_gamma_eq_pi2(), R, 1e-16) == 0);
	free(R);

	{
	double R_correct[9];
	memcpy(R_correct, R_gamma_eq_pi2(), sizeof(R_correct));
	R_mult(R_correct, R_beta_eq_pi2());
	R_mult(R_correct, R_alpha_eq_pi2());
	R = euler_rotation_matrix(M_PI_2, M_PI_2, M_PI_2);
	assert(R_cmp(R_correct, R, 1e-16) == 0);
	free(R);
	}

	R = euler_rotation_matrix(M_PI_4, 0., -M_PI_4);
	assert(R_cmp(R_ident(), R, 1e-16) == 0);
	free(R);

	R = euler_rotation_matrix(M_PI_2, 0., -M_PI_2);
	assert(R_cmp(R_ident(), R, 1e-16) == 0);
	free(R);

	R = euler_rotation_matrix(3. * M_PI_4, 0., -3. * M_PI_4);
	assert(R_cmp(R_ident(), R, 1e-16) == 0);
	free(R);

	R = euler_rotation_matrix(M_PI, 0., -M_PI);
	assert(R_cmp(R_ident(), R, 1e-16) == 0);
	free(R);

	{
	int i;
	for(i = 0; i < 1000000; i++) {
		double alpha = randrange(0., 2. * M_PI);
		double beta = randrange(0., M_PI);
		double gamma = randrange(0., 2. * M_PI);
		double *R1 = euler_rotation_matrix(gamma, beta, alpha);
		R = euler_inv_rotation_matrix(gamma, beta, alpha);
		R_mult(R, R1);
		assert(R_cmp(R_ident(), R, 1e-15) == 0);
		free(R);
		free(R1);
	}
	}

	return 0;
}


static double *vector_from_radec(double ra, double dec)
{
	double *v = malloc(3 * sizeof(*v));

	v[0] = cos(dec) * cos(ra);
	v[1] = cos(dec) * sin(ra);
	v[2] = sin(dec);

	return v;
}


static void radec_from_vector(double *ra, double *dec, double *v)
{
	gsl_vector_view view = gsl_vector_view_array(v, 3);
	vector_direction(&view.vector, dec, ra);
	*dec = M_PI_2 - *dec;
}


static double *vector_rotate(double *v, double *R)
{
	double tmp[] = {0., 0., 0.};
	int i, j;

	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			tmp[i] += R[3 * i + j] * v[j];
	/*fprintf(stderr, "[%8.5g]   [%8.5g %8.5g %8.5g] [%8.5g]\n", tmp[0], R[0], R[1], R[2], v[0]);
	fprintf(stderr, "[%8.5g] = [%8.5g %8.5g %8.5g] [%8.5g]\n", tmp[1], R[3], R[4], R[5], v[1]);
	fprintf(stderr, "[%8.5g]   [%8.5g %8.5g %8.5g] [%8.5g]\n", tmp[2], R[6], R[7], R[8], v[2]);*/
	memcpy(v, tmp, sizeof(tmp));

	return v;
}


static int test_galactic_rotation_matrix(void)
{
	/* does sky_equatorial_to_galactic_rot_matrix() rotate the
	 * equatorial co-ordinates of the galactic co-ordinate system's x
	 * axis to (0, 0) and z axis to (0, pi/2)? */
	{
	/* components of galactic co-ordinate system's x axis in equatorial
	 * co-ordinates */
	double *v = vector_from_radec(SKY_MW_X_J2000_RA_RAD, SKY_MW_X_J2000_DEC_RAD);
	/* rotation matrix to rotate the equatorial co-ordinate system to
	 * the galactic co-ordinate system */
	double *R = sky_equatorial_to_galactic_rot_matrix();
	double ra, dec;
	/* after the rotation the vector should lie on the x axis (because
	 * it is this co-ordinate system's x axis) */
	vector_rotate(v, R);
	radec_from_vector(&ra, &dec, v);
	free(v);
	/*fprintf(stderr, "X rotated to (%g pi, %g pi), expected (0., 0.)\n", ra / M_PI, dec / M_PI);*/
	assert(fabs(ra - 0.) < 1e-3 && fabs(dec - 0.) < 1e-3);
	/* components of the galactic co-ordinate system's z axis in
	 * equatorial co-ordinates */
	v = vector_from_radec(SKY_MW_Z_J2000_RA_RAD, SKY_MW_Z_J2000_DEC_RAD);
	/* after the rotation the vector should lie on the z axis (because
	 * it is this co-ordinate system's z axis) */
	vector_rotate(v, R);
	radec_from_vector(&ra, &dec, v);
	free(v);
	free(R);
	/*fprintf(stderr, "Z rotated to (%g pi, %g pi), expected (0., 0.5 pi)\n", ra / M_PI, dec / M_PI);*/
	/* if dec is correct then ra is degenerate, so we don't check it */
	assert(fabs(dec - M_PI_2) < 1e-14);
	}

	/* are sky_equatorial_to_galactic_rot_matrix() and
	 * sky_galactic_to_equatorial_rot_matrix() each other's inverse? */
	{
	double *R = sky_equatorial_to_galactic_rot_matrix();
	double *Rinv = sky_galactic_to_equatorial_rot_matrix();
	assert(R_cmp(R_ident(), R_mult(R, Rinv), 1e-15) == 0);
	free(R);
	free(Rinv);
	}

	return 0;
}


/*
 * ============================================================================
 *
 *                           Wigner D Matrix Tests
 *
 * ============================================================================
 */


/*
 * do rotations about the z axis produce diagonal D matrixes?
 */

static int wigner_D_test1(void)
{
	struct sh_series *series = sh_series_new(10, 0);
	double omega;
	int l, m, m_prime;

	for(omega = 0.125; omega < 7.0; omega += 0.125) {
		double *R = euler_rotation_matrix(omega, 0., 0.);
		struct sh_series_rotation_plan *plan = sh_series_rotation_plan_new(series, R);
		free(R);

		for(l = 1; l <= (int) series->l_max; l++)
			for(m = -l; m <= +l; m++)
				for(m_prime = -l; m_prime <= +l; m_prime++)
					if(m_prime != m)
						assert(sh_series_rotation_plan_wigner_D(plan, l, m, m_prime) == 0.);

		sh_series_rotation_plan_free(plan);
	}

	sh_series_free(series);

	return 0;
}


/*
 * generate an impulse at some co-ordinates, rotate by some arbitrary
 * angles, compare to correct answer by generating another impulse at the
 * new co-ordinates.  testing with l_max = 1 validates the D1 function in
 * the Wigner D matrix code, testing with l_max > 1 validates the recursion
 * relation in the Wigner D matrix code.
 */

static int wigner_D_test2(unsigned int l_max)
{
	int i;

	for(i = 0; i < 100000; i++) {
		double theta = randrange(0., M_PI);
		double phi = randrange(0., 2. * M_PI);
		double dtheta = randrange(0., M_PI) - theta;
		double dphi = randrange(0., 2. * M_PI) - phi;
		struct sh_series *series = sh_series_impulse(l_max, theta, phi);
		struct sh_series *result = sh_series_new(series->l_max, 0);
		struct sh_series *correct = sh_series_impulse(series->l_max, theta + dtheta, phi + dphi);
		double *R = euler_rotation_matrix(-phi, dtheta, phi + dphi);
		struct sh_series_rotation_plan *plan = sh_series_rotation_plan_new(series, R);

		free(R);
		sh_series_rotate(result, series, plan);
		sh_series_rotation_plan_free(plan);
		sh_series_free(series);

		if(diagnostics_rms_error(result, correct) >= 1e-10) {
			fprintf(stderr, "(theta,phi)=(%.16g,%.16g) (dtheta,dphi)=(%.16g,%.16g)\n", theta, phi, dtheta, dphi);
			sh_series_add(result, -1., correct);
			sh_series_print(stderr, result);
			assert(0);
		}

		sh_series_free(correct);
		sh_series_free(result);
	}

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
	/*
	 * self tests
	 */

	{
	int i;
	for(i = 0; i < 1000000; i++) {
		double ra = randrange(0., 2 * M_PI);
		double dec = randrange(-M_PI_2, +M_PI_2);
		double ra1, dec1;
		double *v = vector_from_radec(ra, dec);
		radec_from_vector(&ra1, &dec1, v);
		free(v);
		if(ra1 < 0.)
			ra1 += 2. * M_PI;
		/*fprintf(stderr, "(%.17g,%.17g) (%.17g,%.17g)\n", ra, dec, ra1, dec1);*/
		assert(fabs(ra - ra1) < 1e-13 && fabs(dec - dec1) < 1e-15);
	}
	}

	{
	int l, m;
	/* should contain only m=0 components */
	struct sh_series *series = sh_series_impulse(10, 0.0, 1.0);
	for(l = 0; l <= (int) series->l_max; l++)
		for(m = -l; m <= l; m++)
			if(m != 0)
				assert(sh_series_get(series, l, m) == 0.);
	sh_series_free(series);
	series = sh_series_impulse(10, M_PI, 2.0);
	for(l = 0; l <= (int) series->l_max; l++)
		for(m = -l; m <= l; m++)
			if(m != 0)
				assert(sh_series_get(series, l, m) == 0.);
	sh_series_free(series);
	}

	/*
	 * test of library
	 */

	assert(test_euler_rotation_matrix() == 0);
	assert(test_galactic_rotation_matrix() == 0);
	assert(wigner_D_test1() == 0);
	assert(wigner_D_test2(1) == 0);
	assert(wigner_D_test2(6) == 0);

	exit(0);
}
