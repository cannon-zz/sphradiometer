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


static struct sh_series *random_sh_series(int l_max, int polar)
{
	struct sh_series *series = sh_series_new(l_max, polar);
	int l, m;

	for(l = 0; l <= l_max; l++)
		for(m = (polar ? 0 : -l); m <= (polar ? 0 : +l); m++)
			sh_series_set(series, l, m, random() / (double) RAND_MAX + I * random() / (double) RAND_MAX);

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


static int sh_series_cmp(const struct sh_series *a, const struct sh_series *b)
{
	/* coefficients must agree up to the smaller of the two, and then
	 * be 0 in the larger of the two */
	int a_l_max = a->l_max;
	int a_m_max = a->polar ? 0 : a_l_max;
	int b_l_max = b->l_max;
	int b_m_max = b->polar ? 0 : b_l_max;
	int l_max = a_l_max < b_l_max ? a_l_max : b_l_max;
	int m_max = a_m_max < b_m_max ? a_m_max : b_m_max;
	int l, m;

	/* check the l for which they both have coefficients */
	for(l = 0; l <= l_max; l++) {
		/* check +ve and -ve m.  wastes cpu cycles for m=0 but who
		 * cares */
		for(m = 0; m <= (l < m_max ? l : m_max); m++)
			if(sh_series_get(a, l, m) != sh_series_get(b, l, m) || sh_series_get(a, l, -m) != sh_series_get(b, l, -m))
				return 1;
		/* any extra elements in a or b must be 0 */
		for(m = m_max + 1; m <= (l < a_m_max ? l : a_m_max); m++)
			if(sh_series_get(a, l, m) != 0.)
				return 1;
		for(m = m_max + 1; m <= (l < b_m_max ? l : b_m_max); m++)
			if(sh_series_get(b, l, m) != 0.)
				return 1;
	}

	/* any extra elements in a or b must be 0 */
	for(l = l_max + 1; l <= a_l_max; l++)
		for(m = 0; m <= (l < a_m_max ? l : a_m_max); m++)
			if(sh_series_get(a, l, m) != 0.)
				return 1;
	for(l = l_max + 1; l <= b_l_max; l++)
		for(m = 0; m <= (l < b_m_max ? l : b_m_max); m++)
			if(sh_series_get(b, l, m) != 0.)
				return 1;

	/* they're equal */
	return 0;
}


static int sh_series_array_cmp(const struct sh_series_array *a, const struct sh_series_array *b)
{
	/* they need to pass sh_series_cmp() up to the shorter of the two
	 * series, and then the longer must be zero beyond that */
	int n = a->n < b->n ? a->n : b->n;
	int i;
	unsigned j;

	for(i = 0; i < n; i++)
		if(sh_series_cmp(&a->series[i], &b->series[i]))
			return 1;

	for(i = n; i < a->n; i++) {
		struct sh_series *series = &a->series[i];
		for(j = 0; j < sh_series_length(series->l_max, series->polar); j++)
			if(series->coeff[j] != 0.)
				return 1;
	}

	for(i = n; i < b->n; i++) {
		struct sh_series *series = &b->series[i];
		for(j = 0; j < sh_series_length(series->l_max, series->polar); j++)
			if(series->coeff[j] != 0.)
				return 1;
	}

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


static complex double Y_10_03(double theta, double phi, void *nul)
{
	return -3.0 / 256 * sqrt(5005 / M_PI) * cexp(I * 3 * phi) * pow(sin(theta), 3) * (323 * pow(cos(theta), 7) - 357 * pow(cos(theta), 5) + 105 * pow(cos(theta), 3) - 7 * cos(theta));
}


static int test_projection(void)
{
	struct sh_series *test = sh_series_new(20, 0);
	struct sh_series *exact = sh_series_new(20, 0);

	sh_series_from_func(test, Y_00_00, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 0, 0, 1);
	fprintf(stderr, "Projection of Y_{0,0} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_from_func(test, Y_07_01, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 7, 1, 1);
	fprintf(stderr, "Projection of Y_{7,1} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_from_func(test, Y_08_n7, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 8, -7, 1);
	fprintf(stderr, "Projection of Y_{8,-7} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

	sh_series_from_func(test, Y_10_03, NULL);
	sh_series_zero(exact);
	sh_series_set(exact, 10, 3, 1);
	fprintf(stderr, "Projection of Y_{10,3} rms error = %g\n", diagnostics_rms_error(test, exact) / (4 * M_PI));
	sh_series_print(stderr, test);

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
	 * test projections of real-valued functions
	 */

	test_projection();

	exit(0);
}
