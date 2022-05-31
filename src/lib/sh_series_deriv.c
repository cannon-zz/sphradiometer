
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
#include <math.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                          Second Order Derivatives
 *
 * ============================================================================
 */


/*
 * Replace an sh_series object with its Laplacian.  Spherical harmonics are
 * eigenfunctions of the Laplacian,
 *
 * 	\grad^{2} Y_{lm}(\hat{s}) = -l (l + 1) Y_{lm}(\hat{s}).
 *
 * Each coefficient of the series is multiplyied by the corresponding
 * eigenvalue.
 */


struct sh_series *sh_series_laplacian(struct sh_series *series)
{
	complex double *c = series->coeff;
	int l = 0, m = 0;

	while(1) {
		*(c++) *= -l * (l + 1);
		if(++l > (int) series->l_max) {
			if(++m > (int) series->l_max)
				m = -(int) series->l_max;
			else if(m == 0)
				break;
			l = abs(m);
		}
	}

	return series;
}


/*
 * Replace an sh_series object with its inverse Laplacian by dividing each
 * coefficient by the corresponding eigenvalue.  The (0,0) coefficient is
 * undefined in this operation (being the arbitrary integration constant),
 * and this funtion sets it to 0.
 */


struct sh_series *sh_series_invlaplacian(struct sh_series *series)
{
	complex double *c = series->coeff;
	int l = 1, m = 0;

	*(c++) = 0.0;
	/* done? */
	if(series->l_max == 0)
		return series;

	while(1) {
		*(c++) /= -l * (l + 1);
		if(++l > (int) series->l_max) {
			if(++m > (int) series->l_max)
				m = -(int) series->l_max;
			else if(m == 0)
				break;
			l = abs(m);
		}
	}

	return series;
}

