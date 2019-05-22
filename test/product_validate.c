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
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                         Random Function on Sphere
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


/*
 * ============================================================================
 *
 *                    Confirm the Product of Two Functions
 *
 * ============================================================================
 */


static int test1(unsigned int lmax1, int polar1, unsigned int lmax2, int polar2)
{
	struct sh_series *series1 = random_sh_series(lmax1, polar1);
	struct sh_series *series2 = random_sh_series(lmax2, polar2);
	struct sh_series *prod = sh_series_new(lmax1 + lmax2, polar1 && polar2);
	struct sh_series_product_plan *plan = NULL;
	double max_err = 0.;
	int i;
	fprintf(stderr, "(l_max_1, polar_1) = (%u, %d); ", lmax1, polar1);
	fprintf(stderr, "(l_max_2, polar_2) = (%u, %d)\n", lmax2, polar2);
	if(!series1 || !series2 || !prod) {
		fprintf(stderr, "%s(): allocation failure\n", __func__);
		goto error;
	}
	plan = sh_series_product_plan_new(prod, series1, series2);
	if(!plan) {
		fprintf(stderr, "%s(): plan creation failure\n", __func__);
		goto error;
	}
	if(!sh_series_product(prod, series1, series2, plan)) {
		fprintf(stderr, "%s(): sh_series_product() failed\n", __func__);
		goto error;
	}

	/* pick a few places on the sphere and confirm that the product
	 * equals the product of the two functions */
	for(i = 0; i < 100; i++) {
		double theta = randrange(0., M_PI);
		double phi = randrange(0., 2. * M_PI);
		complex double val1 = sh_series_eval(series1, theta, phi);
		complex double val2 = sh_series_eval(series2, theta, phi);
		complex double valprod = sh_series_eval(prod, theta, phi);
		double err = cabs(val1 * val2 - valprod);
		if(err > max_err) {
			fprintf(stderr, "%g+i*%g * %g+i*%g = %g+i*%g, abs(error) = %g\n", creal(val1), cimag(val1), creal(val2), cimag(val2), creal(valprod), cimag(valprod), err);
			max_err = err;
		}
	}
	if(max_err > 1e-10)
		goto error;

	sh_series_free(series1);
	sh_series_free(series2);
	sh_series_free(prod);
	sh_series_product_plan_free(plan);

	return 0;

error:
	sh_series_free(series1);
	sh_series_free(series2);
	sh_series_free(prod);
	sh_series_product_plan_free(plan);
	return -1;
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
	 * test of library
	 */

	assert(test1(5, 1, 10, 1) == 0);
	assert(test1(5, 0, 10, 0) == 0);
	assert(test1(70, 1, 70, 1) == 0);

	exit(0);
}
