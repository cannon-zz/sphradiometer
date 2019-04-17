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
		for(m = polar ? 0 : -l; m <= polar ? 0 : l; m++)
			sh_series_set(series, l, m, random() / (double) RAND_MAX + I * random() / (double) RAND_MAX);

	return series;
}


static struct sh_series_array *random_sh_series_array(int l_max, int polar, int n)
{
	struct sh_series_array *array = sh_series_array_new(l_max, polar, n);
	int i;

	for(i = 0; i < n; i++) {
		struct sh_series *series = random_sh_series(l_max, polar);
		sh_series_assign(&array->series[i], series);
		sh_series_free(series);
	}

	return array;
}


static int sh_series_cmp(const struct sh_series *a, const struct sh_series *b)
{
	/* coefficients must agree up to the smaller of the two, and then
	 * be 0 in the larger of the two */
	size_t a_len = sh_series_length(a->l_max, a->polar);
	size_t b_len = sh_series_length(b->l_max , b->polar);
	size_t len = a_len < b_len ? a_len : b_len;
	size_t i;

	if(memcmp(a->coeff, b->coeff, len * sizeof(*a->coeff)))
		return 1;
	for(i = len; i < a_len; i++)
		if(a->coeff[i] != 0.)
			return 1;
	for(i = len; i < b_len; i++)
		if(b->coeff[i] != 0.)
			return 1;

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
	 * confirm that sh_series_set_polar() and
	 * sh_series_array_set_polar() preserve the value
	 */

	{
	struct sh_series *a = random_sh_series(57, 1);
	struct sh_series *b = sh_series_copy(a);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_set_polar(a, 0);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_set_polar(a, 1);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	{
	struct sh_series_array *a = random_sh_series_array(57, 1, 15);
	struct sh_series_array *b = sh_series_array_copy(a);
	assert(sh_series_array_cmp(a, b) == 0);
	sh_series_array_set_polar(a, 0);
	assert(sh_series_array_cmp(a, b) == 0);
	sh_series_array_set_polar(a, 1);
	assert(sh_series_array_cmp(a, b) == 0);
	sh_series_array_free(a);
	sh_series_array_free(b);
	}

	/*
	 * confirm that sh_series_assig() preserves the value
	 */

	{
	struct sh_series *a = random_sh_series(57, 1);
	struct sh_series *b = sh_series_new(57, 0);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_assign(b, a);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/*
	 * test projections of real-valued functions
	 */

	test_projection();

	exit(0);
}
