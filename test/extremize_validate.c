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
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <sphradiometer/diagnostics.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                   Tests
 *
 * ============================================================================
 */


static int test_real_maximum(void)
{
	int l_max = 22;

	int i;
	for(i = 0; i < 20; i++) {
		double correct_theta = randrange(0., M_PI);
		double correct_phi = randrange(0., 2. * M_PI);
		struct sh_series *series = sh_series_impulse(l_max, correct_theta, correct_phi);

		double theta, phi;
		double max;

		fprintf(stderr, "expected max %.16g @ theta, phi = %.16g %.16g\n", cabs(sh_series_eval(series, correct_theta, correct_phi)), correct_theta, correct_phi);

		sh_series_neg(series);
		max = -sh_series_real_minimum(series, &theta, &phi);
		sh_series_free(series);

		fprintf(stderr, "estimated theta, phi = %.16g %.16g\n", theta, phi);
		fprintf(stderr, "value at theta, phi = %.16g\n", max);
		fprintf(stderr, "\n");
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
	srandom(time(NULL));

	/*
	 * sh_series_real_maximum()
	 */

	assert(test_real_maximum() == 0);

	exit(0);
}
