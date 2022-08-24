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
#include <unistd.h>

#include <sphradiometer/diagnostics.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/sphradiometer_config.h> /* needed for HAVE_GSL_2_0 */


/*
 * ============================================================================
 *
 *                                  Helpers
 *
 * ============================================================================
 */


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
 *                                Command Line
 *
 * ============================================================================
 */


struct options {
	int benchmark;
};


static struct options parse_command_line(int argc, char * const argv[])
{
	struct options options = {
		.benchmark = 0,
	};
	char opt;

	while((opt = getopt(argc, argv, "bh")) != -1) {
		switch(opt) {
		case 'b':
			options.benchmark = 1;
			break;

		case 'h':
			fprintf(stderr, "usage:\n\t%s [options]\n\nOptions:\n\t-b\tRun in benchmark mode (default = run in test suite mode).\n\t-h\tDisplay this message.\n\n", argv[0]);
			exit(0);

		case '?':
		case ':':
		default:
			fprintf(stderr, "command line error\n");
			exit(1);
		}
	}

	return options;
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


static complex double Y_10_n10(double theta, double phi, void *nul)
{
	return 1. / 1024 * sqrt(969969 / M_PI) * cexp(-I * 10 * phi) * pow(sin(theta), 10);
}


static complex double Y_10_03(double theta, double phi, void *nul)
{
	return -3.0 / 256 * sqrt(5005 / M_PI) * cexp(I * 3 * phi) * pow(sin(theta), 3) * (323 * pow(cos(theta), 7) - 357 * pow(cos(theta), 5) + 105 * pow(cos(theta), 3) - 7 * cos(theta));
}


static complex double Y_10_10(double theta, double phi, void *nul)
{
	return 1. / 1024 * sqrt(969969 / M_PI) * cexp(I * 10 * phi) * pow(sin(theta), 10);
}


static double reY_wrapper(double theta, double phi, void *f)
{
	complex double (*Y)(double, double, void *) = f;

	return creal(Y(theta, phi, NULL));
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
		{10, -10, Y_10_n10},
		{10, 3, Y_10_03},
		{10, +10, Y_10_10},
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


static int test_conj(void)
{
	int i;
	for(i = 0; i < 100; i++) {
		struct sh_series *a = random_sh_series(floor(randrange(0, 21)), random() & 1);
		struct sh_series *b = sh_series_conj(sh_series_copy(a));
		int j;
		assert(b != NULL);
		for(j = 0; j < 100; j++) {
			double theta = randrange(0., M_PI);
			double phi = randrange(0., 2. * M_PI);
			complex double a_val_conj = conj(sh_series_eval(a, theta, phi));
			complex double b_val = sh_series_eval(b, theta, phi);
			double err = cabs(b_val - a_val_conj);
			if(err > 1e-13) {
				fprintf(stderr, "%.17g+I%.17g != %.17g+I%.17g (cabs(error) = %g)\n", creal(b_val), cimag(b_val), creal(a_val_conj), cimag(a_val_conj), err);
				return 1;
			}
		}
		sh_series_free(a);
		sh_series_free(b);
	}

	return 0;
}


static int test_realimag(void)
{
	int i;

	for(i = 0; i < 100; i++) {
		struct sh_series *orig = random_sh_series(floor(randrange(0, 21)), random() & 1);
		struct sh_series *real = sh_series_real(sh_series_copy(orig));
		struct sh_series *imag = sh_series_imag(sh_series_copy(orig));
		struct sh_series *sum;
		int j;
		double err;
		assert(real != NULL);
		assert(imag != NULL);
		sum = sh_series_copy(real);
		assert(sum != NULL);
		/* does s = real(s) + I * imag(s) ? */
		sh_series_add(sum, I*1., imag);
		err = diagnostics_rms_error(orig, sum);
		if(err > 1e-15) {
			fprintf(stderr, "for (l,polar) = (%u,%d), testing s=real(s)+I*imag(s) rms error = %g\norig:", orig->l_max, orig->polar, err);
			sh_series_print(stderr, orig);fprintf(stderr, "\nreal:");
			sh_series_print(stderr, real);fprintf(stderr, "\nimag:");
			sh_series_print(stderr, imag);fprintf(stderr, "\nsum:");
			sh_series_print(stderr, sum);fprintf(stderr, "\n");
			return 1;
		}
		/* are real(s) and imag(s) correct?  check some random
		 * co-ordinates */
		for(j = 0; j < 100; j++) {
			double theta = randrange(0., M_PI);
			double phi = randrange(0., 2. * M_PI);
			complex double y = sh_series_eval(orig, theta, phi);
			complex double x = sh_series_eval(real, theta, phi);
			if(cabs(x - creal(y)) > 1e-13) {
				fprintf(stderr, "(l,polar)=(%u,%d) @ (theta,phi)=(%g,%g) expected %.17g+I%.17g got %.17g+I*%.17g\nseries was:", orig->l_max, orig->polar, theta, phi, creal(y), 0., creal(x), cimag(x));
				sh_series_print(stderr, orig);fprintf(stderr, "\nreal:");
				sh_series_print(stderr, real);fprintf(stderr, "\nimag:");
				sh_series_print(stderr, imag);fprintf(stderr, "\nsum:");
				sh_series_print(stderr, sum);fprintf(stderr, "\n");
				return 1;
			}
			x = sh_series_eval(imag, theta, phi);
			if(cabs(x - cimag(y)) > 1e-13) {
				fprintf(stderr, "(l,polar)=(%u,%d) @ (theta,phi)=(%g,%g) expected %.17g+I%.17g got %.17g+I*%.17g\nseries was:", orig->l_max, orig->polar, theta, phi, cimag(y), 0., creal(x), cimag(x));
				sh_series_print(stderr, orig);fprintf(stderr, "\nreal:");
				sh_series_print(stderr, real);fprintf(stderr, "\nimag:");
				sh_series_print(stderr, imag);fprintf(stderr, "\nsum:");
				sh_series_print(stderr, sum);fprintf(stderr, "\n");
				return 1;
			}
		}
		sh_series_free(orig);
		sh_series_free(real);
		sh_series_free(imag);
		sh_series_free(sum);
	}

	return 0;
}


static int test_add(void)
{
	int i;
	for(i = 0; i < 10000; i++) {
		struct sh_series *a = random_sh_series(floor(randrange(0, 21)), random() & 1);
		struct sh_series *b = random_sh_series(floor(randrange(0, a->l_max)), (random() & 1) | a->polar);
		complex double z = randrange(1., 10.);
		struct sh_series *c = sh_series_add(sh_series_copy(a), z, b);
		int j;

		for(j = 0; j < 10; j++) {
			double theta = randrange(0., M_PI);
			double phi = randrange(0., 2. * M_PI);
			complex double A_plus_z_times_B = sh_series_eval(a, theta, phi) + z * sh_series_eval(b, theta, phi);
			complex double C = sh_series_eval(c, theta, phi);

			double err = cabs(C - A_plus_z_times_B);
			if(err > 1e-13) {
				fprintf(stderr, "sh_series_add failed:  a: l_max=%d,polar=%d,  b: l_max=%d,polar=%d.\n", a->l_max, a->polar, b->l_max, b->polar);
				return -1;
			}
		}

		sh_series_free(a);
		sh_series_free(b);
		sh_series_free(c);
	}

	return 0;
}


static int test_projection1(void)
{
	struct sh_series *test = sh_series_new(10, 0);
	struct sh_series *realtest = sh_series_new(10, 0);
	struct sh_series *exact = sh_series_new(10, 0);
	struct sh_series *realexact = sh_series_new(10, 0);
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
		{10, -10, Y_10_n10},
		{10, 3, Y_10_03},
		{10, +10, Y_10_10},
		{0, 0, NULL}
	};
	int i, j;

	for(i = 0; tests[i].func; i++) {
		double err;
		sh_series_from_func(test, tests[i].func, NULL);
		sh_series_from_realfunc(realtest, reY_wrapper, tests[i].func);
		sh_series_zero(exact);
		sh_series_zero(realexact);
		sh_series_set(exact, tests[i].l, tests[i].m, 1);
		if(tests[i].l) {
			sh_series_set(realexact, tests[i].l, tests[i].m, 0.5);
			sh_series_set(realexact, tests[i].l, -tests[i].m, (tests[i].m & 1) ? -0.5 : 0.5);
		} else {
			sh_series_set(realexact, tests[i].l, tests[i].m, 1);
		}
		err = diagnostics_rms_error(test, exact);
		if(err > 1e-14) {
			fprintf(stderr, "Projection of Y_{%d,%d} rms error = %g\n", tests[i].l, tests[i].m, err);
			sh_series_clip(test, 1e-12);
			sh_series_print(stderr, test);
			return -1;
		}
		err = diagnostics_rms_error(realtest, realexact);
		if(err > 1e-14) {
			fprintf(stderr, "Projection of Re Y_{%d,%d} rms error = %g\n", tests[i].l, tests[i].m, err);
			sh_series_clip(realtest, 1e-12);
			sh_series_print(stderr, realtest);
			return -1;
		}
		for(j = 0; j < 100; j++) {
			double theta = randrange(0., M_PI);
			double phi = randrange(0., 2. * M_PI);
			err = cabs(sh_series_eval(exact, theta, phi) - tests[i].func(theta, phi, NULL));
			if(err > 5e-10) {
				fprintf(stderr, "sh_series_eval(..., %g, %g) failed for (l,m)=(%d,%d), |err|=%g\n", theta, phi, tests[i].l, tests[i].m, err);
				return -1;
			}
			err = cabs(sh_series_eval(realexact, theta, phi) - sh_series_eval(realtest, theta, phi));
			if(err > 1e-13) {
				fprintf(stderr, "sh_series_from_realfunc() failed for (l,m)=(%d,%d), |err|=%g\n", tests[i].l, tests[i].m, err);
				return -1;
			}
		}
	}

	sh_series_free(test);
	sh_series_free(exact);

	return 0;
}


static int test_projection2(void)
{
	double err;
	unsigned int lmax;

	for(lmax = 1; lmax < 50; lmax++) {
		struct sh_series *orig = random_sh_series(lmax, 0);
		complex double *mesh = sh_series_to_mesh(orig);
		struct sh_series *final = sh_series_new(lmax, 0);
		sh_series_from_mesh(final, mesh);
		err = diagnostics_rms_error(orig, final) / lmax / lmax;
		free(mesh);
		sh_series_free(orig);
		sh_series_free(final);
		if(err > 1e-12) {
			fprintf(stderr, "sh_series_from_mesh(sh_series_to_mesh(x)) != x for lmax=%u. |err| = %g\n", lmax, err);
			return -1;
		}
	}
	return 0;
}


#ifdef HAVE_GSL_2_0
static complex double interp_wrapper(double theta, double phi, void *interp)
{
	return sh_series_eval_interp((struct sh_series_eval_interp *) interp, theta, phi);
}


static int test_interpolation(void)
{
	unsigned int lmax;

	for(lmax = 1; lmax < 50; lmax++) {
		/* construct a random sh_series and build an interpolator
		 * for its pixel domain representation */
		struct sh_series *orig = random_sh_series(lmax, 0);
		struct sh_series_eval_interp *interp = sh_series_eval_interp_new(orig);
		struct sh_series *final = sh_series_new(2 * lmax, 0);
		double err;

		/* check a few properties of the interpolator.  can we
		 * evaluate it at the poles, is it periodic in phi, etc. */

		double complex x, y;
		(void) sh_series_eval_interp(interp, 0., 1.);
		(void) sh_series_eval_interp(interp, M_PI, 1.);
		x = sh_series_eval_interp(interp, M_PI / 2., 0.);
		y = sh_series_eval_interp(interp, M_PI / 2., 2. * M_PI);
		if(x != y) {
			fprintf(stderr, "for lmax=%d theta=pi/2 @phi=0 x=%.16g+i*%.16g @phi=2pi x=%.16g+i*%.16g\n", lmax, creal(x), cimag(x), creal(y), cimag(y));
			return -1;
		}
		y = sh_series_eval_interp(interp, M_PI / 2., -4. * M_PI);
		if(x != y) {
			fprintf(stderr, "for lmax=%d theta=pi/2 @phi=0 x=%.16g+i*%.16g @phi=-4pi x=%.16g+i*%.16g\n", lmax, creal(x), cimag(x), creal(y), cimag(y));
			return -1;
		}

		sh_series_from_func(final, interp_wrapper, interp);
		sh_series_eval_interp_free(interp);
		err = diagnostics_rms_error(orig, final) / lmax;
		/* FIXME:  this is surprisingly large, but it seems to be
		 * correct ... investigate */
		if(err > 0.03) {
			fprintf(stderr, "orig:\n");
			sh_series_print(stderr, orig);
			fprintf(stderr, "final:\n");
			sh_series_print(stderr, final);
			fprintf(stderr, "sh_series_from_func(sh_series_eval_interp(x)) != x for lmax=%u. |err| = %g\n", lmax, err);
			return -1;
		}
		sh_series_free(final);
		sh_series_free(orig);
	}
	return 0;
}
#endif


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	struct options options = parse_command_line(argc, argv);

	srandom(time(NULL));

	if(options.benchmark) {
		for(int i = 0; i < 1000; i++) {
			assert(test_projection1() == 0);
			assert(test_projection2() == 0);
		}
		exit(0);
	}

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

	/* non overlapping memory */
	{
	int i;
	for(i = 0; i < 100000; i++) {
		struct sh_series *dst = random_sh_series(random() & 0x7f, random() & 0x1);
		struct sh_series *src = sh_series_copy(dst);
		sh_series_resize(src, random() & 0x7f);
		if(random() & 0x1)
			sh_series_set_polar(src, !src->polar);
		/* src and dst now have random, possibly different, sizes,
		 * but dst is guaranteed to be able to accept all of the
		 * non-zero coefficients from src, so they should be
		 * functionally identical after an assign operation */
		assert(sh_series_assign(dst, src) != NULL);
		assert(sh_series_cmp(dst, src) == 0);
		sh_series_free(dst);
		sh_series_free(src);
	}
	}

	/* overlapping memory */
	{
	int i;
	for(i = 0; i < 100000; i++) {
		struct sh_series *dst = random_sh_series(random() & 0x7f, random() & 0x1);
		struct sh_series *src = sh_series_copy(dst);
		struct sh_series wrapper;
		sh_series_resize(src, random() & 0x7f);
		if(random() & 0x1)
			sh_series_set_polar(src, !src->polar);
		/* src and dst now have random, possibly different, sizes,
		 * but dst is guaranteed to be able to accept all of the
		 * non-zero coefficients from src, so they should be
		 * functionally identical after an assign operation */
		/* we're now going to copy src's contents into the memory
		 * used to hold dst's coefficients, then assign from the
		 * former to the latter */
		if(sh_series_length(src->l_max, src->polar) > sh_series_length(dst->l_max, dst->polar)) {
			/* make room in dst */
			dst->coeff = realloc(dst->coeff, sh_series_length(src->l_max, src->polar) * sizeof(*dst->coeff));
			assert(dst->coeff != NULL);
		}
		memcpy(dst->coeff, src->coeff, sh_series_length(src->l_max, src->polar) * sizeof(*dst->coeff));
		wrapper = *src;
		wrapper.coeff = dst->coeff;
		assert(sh_series_assign(dst, &wrapper) != NULL);
		assert(sh_series_cmp(dst, src) == 0);
		sh_series_free(dst);
		sh_series_free(src);
	}
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
	assert(test_conj() == 0);
	assert(test_realimag() == 0);
	assert(test_add() == 0);
	assert(test_projection1() == 0);
	assert(test_projection2() == 0);
#ifdef HAVE_GSL_2_0
	assert(test_interpolation() == 0);
#endif

	exit(0);
}
