/*
 * Copyright (C) 2021  Kipp C. Cannon
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


#include <stdio.h>


#include <gsl/gsl_sf_legendre.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


/*
 * GSL's (l,m) ordering is
 *
 * (0, 0),
 * (1, 0), (1, 1),
 * (2, 0), (2, 1), (2, 2)
 * (3, 0), (3, 1), (3, 2), (3, 3)
 * ...
 */


int main(int argc, char *argv[])
{
	int l_max = 10;
	int l, m;

	/* confirm that GSL's coefficient ordering is as we believe it to
	 * be */

	for(l = 0; l <= l_max; l++)
		for(m = -l; m <= l; m++)
			if((int) gsl_sf_legendre_array_index(l, m) != l * (l + 1) / 2 + m) {
				fprintf(stderr, "gsl_sf_legendre_array_index(%d, %d) = %zu but expected %d\n", l, m, gsl_sf_legendre_array_index(l, m), l * (l + 1) / 2 + m);
				exit(1);
			}

	return 0;
}
