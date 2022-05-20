/*
 * Copyright (C) 2006--2009,2019  Kipp C. Cannon
 * Copyright (C) 2019  Takuya Tsutsui
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


#include <complex.h>
#include <stdio.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                  File I/O
 *
 * ============================================================================
 */


/*
 * Print the non-zero coefficients in an sh_series object.
 *
 * This function is for diagnostic purposes.  It is not intended to be used
 * for production I/O applications.  The inverse, a parsing function, is
 * avaiable in the Python library.
 */


void sh_series_print(FILE *f, const struct sh_series *series)
{
	int l, m;
	int lines = 0;

	for(l = 0; l <= (int) series->l_max; l++)
		for(m = series->polar ? 0 : -l; m <= (series->polar ? 0 : l); m++) {
			complex double c = sh_series_get(series, l, m);
			if(c != 0.0) {
				fprintf(f, "(%d,%d) = %.17g + I %.17g\n", l, m, creal(c), cimag(c));
				lines++;
			}
		}
	/* if all coefficients are zero, print at least something to
	 * simplify code that tries to parse this output */
	if(!lines)
		fprintf(f, "(%d,%d) = %.17g + I %.17g\n", 0, 0, 0.0, 0.0);
}
