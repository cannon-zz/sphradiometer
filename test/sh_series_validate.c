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
#include <time.h>

#include <sphradiometer/diagnostics.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                  Helpers
 *
 * ============================================================================
 */


static double randrange(double lo, double hi)
{
	return lo + (double) random() * (hi - lo) / RAND_MAX;
}


static struct sh_series *random_sh_series(int l_max, int polar)
{
	struct sh_series *series = sh_series_new(l_max, polar);
	int l, m;

	for(l = 0; l <= l_max; l++)
		for(m = (polar ? 0 : -l); m <= (polar ? 0 : +l); m++)
			sh_series_set(series, l, m, randrange(0., 1.) + I * randrange(0., 1.));

	return series;
}


static struct sh_series_array *random_sh_series_array(int n, int l_max, int polar)
{
	struct sh_series_array *array = sh_series_array_new(n, l_max, polar);
	int i;

	for(i = 0; i < n; i++) {
		struct sh_series *series = random_sh_series(l_max, polar);
		assert(sh_series_assign(&array->series[i], series) != NULL);
		sh_series_free(series);
	}

	return array;
}


static int sh_series_array_cmp(const struct sh_series_array *a, const struct sh_series_array *b)
{
	/* they need to pass sh_series_cmp() up to the shorter of the two
	 * series, and then the longer must be zero beyond that */
	int n = a->n < b->n ? a->n : b->n;
	int i;

	for(i = 0; i < n; i++)
		if(sh_series_cmp(&a->series[i], &b->series[i]))
			return 1;

	for(i = n * a->stride; i < a->n * a->stride; i++)
		if(a->coeff[i] != 0.)
			return 1;

	for(i = n * b->stride; i < b->n * b->stride; i++)
		if(b->coeff[i] != 0.)
			return 1;

	return 0;
}


/*
 * ============================================================================
 *
 *                              sh_series Tests
 *
 * ============================================================================
 */


/*
 * Test projection.  The Y_{lm} test functions have been copied from
 * Wikipedia.
 */


static complex double Y_00_00(double theta, double phi, void *nul)
{
	return 1.0 / 2 * sqrt(1 / M_PI);
}


static complex double Y_07_01(double theta, double phi, void *nul)
{
	return -1.0 / 64 * sqrt(105 / (2 * M_PI)) * cexp(I * phi) * sin(theta) * (429 * pow(cos(theta), 6) - 495 * pow(cos(theta), 4) + 135 * pow(cos(theta), 2) - 5);
}


static complex double Y_08_n7(double theta, double phi, void *nul)
{
	return 3.0 / 64 * sqrt(12155 / (2 * M_PI)) * cexp(I * -7 * phi) * pow(sin(theta), 7) * cos(theta);
}


static complex double Y_08_n6(double theta, double phi, void *nul)
{
	return 1.0 / 128 * sqrt(7293 / M_PI) * cexp(I * -6 * phi) * pow(sin(theta), 6) * (15 * pow(cos(theta), 2) - 1);
}


static complex double Y_08_06(double theta, double phi, void *nul)
{
	return 1.0 / 128 * sqrt(7293 / M_PI) * cexp(I * 6 * phi) * pow(sin(theta), 6) * (15 * pow(cos(theta), 2) - 1);
}


static complex double Y_08_07(double theta, double phi, void *nul)
{
	return -3.0 / 64 * sqrt(12155 / (2 * M_PI)) * cexp(I * 7 * phi) * pow(sin(theta), 7) * cos(theta);
}


static complex double Y_10_03(double theta, double phi, void *nul)
{
	return -3.0 / 256 * sqrt(5005 / M_PI) * cexp(I * 3 * phi) * pow(sin(theta), 3) * (323 * pow(cos(theta), 7) - 357 * pow(cos(theta), 5) + 105 * pow(cos(theta), 3) - 7 * cos(theta));
}


static int test_evaluation1(void)
{
	struct {
		int l;
		int m;
		complex double (*func)(double, double, void *);
	} tests[] = {
		{0, 0, Y_00_00},
		{7, 1, Y_07_01},
		{8, -7, Y_08_n7},
		{8, -6, Y_08_n6},
		{8, +6, Y_08_06},
		{8, +7, Y_08_07},
		{10, 3, Y_10_03},
		{0, 0, NULL}
	};
	int i, j;

	for(i = 0; tests[i].func; i++)
		for(j = 0; j < 100000; j++) {
			double theta = randrange(0., M_PI);
			double phi = randrange(0., 2. * M_PI);
			complex double test = sh_series_Y(tests[i].l, tests[i].m, theta, phi);
			complex double testconj = sh_series_Yconj(tests[i].l, tests[i].m, theta, phi);
			complex double exact = tests[i].func(theta, phi, NULL);
			double err = cabs(test - exact);
			if(err > 1e-9) {
				fprintf(stderr, "Y_{%d,%d}(%g, %g) failure:  expected %.16g+i*%.16g, got %.16g+i*%.16g, |err| = %g\n", tests[i].l, tests[i].m, theta, phi, creal(exact), cimag(exact), creal(test), cimag(test), err);
				return -1;
			}
			if(cabs(test - conj(testconj)) > 1e-17) {
				fprintf(stderr, "Y^*_{%d,%d}(%g, %g) != Yconj_{%d,%d}(%g, %g), got expected %.16g+i*%.16g, got %.16g+i*%.16g, |err| = %g\n", tests[i].l, tests[i].m, theta, phi, tests[i].l, tests[i].m, theta, phi, creal(exact), cimag(exact), creal(test), cimag(test), err);
				return -1;
			}
		}

	return 0;
}


static int test_evaluation2(void)
{
	int i;

	for(i = 0; i < 100000; i++) {
		int lmax = floor(randrange(0, 21));
		int m = floor(randrange(-lmax, +lmax + 1));
		complex double array[lmax - abs(m) + 1];
		complex double arrayconj[lmax - abs(m) + 1];
		double theta = randrange(0., M_PI);
		double phi = randrange(0., 2. * M_PI);
		int l;

		sh_series_Y_array(array, lmax, m, theta, phi);
		sh_series_Yconj_array(arrayconj, lmax, m, theta, phi);

		for(l = abs(m); l <= lmax; l++) {
			complex double expected = sh_series_Y(l, m, theta, phi);
			complex double got = array[l - abs(m)];
			double err = cabs(expected - got);
			if(err > 1e-13) {
				fprintf(stderr, "array Y_{%d,%d}(%g, %g) for lmax=%d failure:  expected %.16g+i*%.16g, got %.16g+i*%.16g, |err| = %g\n", l, m, theta, phi, lmax, creal(expected), cimag(expected), creal(got), cimag(got), err);
				return -1;
			}
			if(cabs(array[l - abs(m)] - conj(arrayconj[l - abs(m)])) > 1e-16) {
				fprintf(stderr, "array Yconj_{%d,%d}(%g, %g) for lmax=%d failed\n", l, m, theta, phi, lmax);
				return -1;
			}
		}
	}

	return 0;
}


static int test_projection(void)
{
	struct sh_series *test = sh_series_new(20, 0);
	struct sh_series *exact = sh_series_new(20, 0);
	struct {
		int l;
		int m;
		complex double (*func)(double, double, void *);
	} tests[] = {
		{0, 0, Y_00_00},
		{7, 1, Y_07_01},
		{8, -7, Y_08_n7},
		{8, -6, Y_08_n6},
		{8, +6, Y_08_06},
		{8, +7, Y_08_07},
		{10, 3, Y_10_03},
		{0, 0, NULL}
	};
	int i, j;

	for(i = 0; tests[i].func; i++) {
		double err;
		sh_series_from_func(test, tests[i].func, NULL);
		sh_series_zero(exact);
		sh_series_set(exact, tests[i].l, tests[i].m, 1);
		err = diagnostics_rms_error(test, exact);
		if(err > 1e-14) {
			fprintf(stderr, "Projection of Y_{%d,%d} rms error = %g\n", tests[i].l, tests[i].m, err);
			sh_series_print(stderr, test);
			return -1;
		}
		for(j = 0; j < 100; j++) {
			double theta = randrange(0., M_PI);
			double phi = randrange(0., 2. * M_PI);
			err = cabs(sh_series_eval(exact, theta, phi) - tests[i].func(theta, phi, NULL));
			if(err > 1e-13) {
				fprintf(stderr, "sh_series_eval(..., %g, %g) failed for (l,m)=(%d,%d), |err|=%g\n", theta, phi, tests[i].l, tests[i].m, err);
				return -1;
			}
		}
	}

	sh_series_free(test);
	sh_series_free(exact);

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
	srandom(time(NULL));

	/*
	 * confirm that sh_series_resize() preserves the value (of
	 * coefficients that survive)
	 */

	/* source not azimuthally symmetric */
	/* first test growing l_max */
	{
	struct sh_series *a = random_sh_series(8, 0);
	struct sh_series *b = sh_series_copy(a);
	assert(sh_series_cmp(a, b) == 0);
	a = sh_series_resize(a, 15);
	assert(a != NULL);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}
	/* now rely on that to make a series with 0'ed coefficients for use
	 * in testing shrinking l_max */
	{
	struct sh_series *a = random_sh_series(8, 0);
	a = sh_series_resize(a, 15);
	struct sh_series *b = sh_series_copy(a);
	assert(sh_series_cmp(a, b) == 0);
	a = sh_series_resize(a, 8);
	assert(a != NULL);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/* azimuthally symmetric */
	/* first test growing l_max */
	{
	struct sh_series *a = random_sh_series(8, 1);
	struct sh_series *b = sh_series_copy(a);
	assert(sh_series_cmp(a, b) == 0);
	a = sh_series_resize(a, 15);
	assert(a != NULL);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}
	/* now rely on that to make a series with 0'ed coefficients for use
	 * in testing shrinking l_max */
	{
	struct sh_series *a = random_sh_series(8, 1);
	a = sh_series_resize(a, 15);
	struct sh_series *b = sh_series_copy(a);
	assert(sh_series_cmp(a, b) == 0);
	a = sh_series_resize(a, 8);
	assert(a != NULL);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/*
	 * confirm that sh_series_assign() preserves the value
	 */

	/* source and destination both not azimuthally symmetric */
	{
	struct sh_series *a = random_sh_series(57, 0);
	struct sh_series *b = sh_series_new(57, 0);
	sh_series_assign(b, a);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/* source azimuthally symmetric, destination not */
	{
	struct sh_series *a = random_sh_series(57, 1);
	struct sh_series *b = sh_series_new(57, 0);
	sh_series_assign(b, a);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/* test with overlapping memory, both azimuthally symmetric */
	{
	struct sh_series *a = random_sh_series(7, 1);
	struct sh_series *b = sh_series_copy(a);
	struct sh_series c = *a;
	c.coeff = a->coeff = realloc(a->coeff, 2 * sh_series_length(c.l_max, c.polar) * sizeof(*c.coeff));
	assert(c.coeff != NULL);
	c.coeff += 4;
	sh_series_assign(&c, a);
	assert(sh_series_cmp(b, &c) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/* test with overlapping memory, src azimuthal, dst not */
	{
	struct sh_series *a = random_sh_series(15, 1);
	struct sh_series *b = sh_series_copy(a);
	struct sh_series c = *a;
	c.polar = 0;
	c.coeff = a->coeff = realloc(a->coeff, 2 * sh_series_length(c.l_max, c.polar) * sizeof(*c.coeff));
	assert(c.coeff != NULL);
	c.coeff += 6;
	sh_series_assign(&c, a);
	assert(sh_series_cmp(b, &c) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/*
	 * confirm that sh_series_set_polar() and
	 * sh_series_array_set_polar() preserve the value
	 */

	{
	struct sh_series *a = random_sh_series(7, 1);
	struct sh_series *b = sh_series_copy(a);
	assert(sh_series_cmp(a, b) == 0);
	assert(sh_series_set_polar(a, 0) != NULL);
	assert(sh_series_cmp(a, b) == 0);
	assert(sh_series_set_polar(a, 1) != NULL);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	{
	struct sh_series_array *a = random_sh_series_array(15, 7, 1);
	struct sh_series_array *b = sh_series_array_copy(a);
	assert(sh_series_array_cmp(a, b) == 0);
	assert(sh_series_array_set_polar(a, 0) != NULL);
	assert(sh_series_array_cmp(a, b) == 0);
	assert(sh_series_array_set_polar(a, 1) != NULL);
	assert(sh_series_array_cmp(a, b) == 0);
	sh_series_array_free(a);
	sh_series_array_free(b);
	}

	/*
	 * other test functions
	 */

	assert(test_evaluation1() == 0);
	assert(test_evaluation2() == 0);
	assert(test_projection() == 0);

	exit(0);
}
