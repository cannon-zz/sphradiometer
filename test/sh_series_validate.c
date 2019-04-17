#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
	 * confirm that sh_series_set_polar() preserves the value
	 */

	{
	struct sh_series *a = random_sh_series(57, 1);
	struct sh_series *b = sh_series_copy(a);
	sh_series_set_polar(a, 0);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_set_polar(a, 1);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	/*
	 * confirm that sh_series_assig() preserves the value
	 */

	{
	struct sh_series *a = random_sh_series(57, 1);
	struct sh_series *b = sh_series_new(57, 0);
	sh_series_assign(b, a);
	assert(sh_series_cmp(a, b) == 0);
	sh_series_free(a);
	sh_series_free(b);
	}

	exit(0);
}
