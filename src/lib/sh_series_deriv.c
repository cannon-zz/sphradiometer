
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
 * Replace an sh_series object with the Laplacian of the function it
 * describes.  Spherical harmonics are eigenfunctions of the Laplacian,
 *
 * 	\grad^{2} Y_{lm}(\hat{s}) = -l (l + 1) Y_{lm}(\hat{s}).
 *
 * The Laplacian is computed by multiplying each coefficient of the input
 * by the corresponding integer eigenvalue.  This is a fast operation.
 *
 * The explicit form of the operator is
 *
 *     1     [ \partial       (           \partial       ) ]
 * --------- [ -------------- ( sin theta -------------- ) ] +
 * sin theta [ \partial theta (           \partial theta ) ]
 *
 *		     1      \partial^2
 *		----------- ---------------
 *		sin^2 theta \partial \phi^2
 *
 * NOTE that because the spherical harmonic basis functions are
 * eigenfunctions of this operator, the 1/sin(theta) factors that appear in
 * this expression do not cause difficulties the way they do when
 * evaluating the divergence or gradient.  In effect, they are guaranteed
 * to cancel terms in the numerator at the poles, and this cancellation
 * occurs term-by-term for each spherical harmonic basis function
 * individually.
 */


struct sh_series *sh_series_laplacian(struct sh_series *series)
{
	complex double *c = series->coeff;
	int l = 0, m = 0;

	while(1) {
		*(c++) *= -l * (l + 1);
		if(++l > (int) series->l_max) {
			if(series->polar)
				break;
			else if(++m > (int) series->l_max)
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

	/* l = 0 */
	*(c++) = 0.0;

	/* done? */
	if(series->l_max == 0)
		return series;

	while(1) {
		*(c++) /= -l * (l + 1);
		if(++l > (int) series->l_max) {
			if(series->polar)
				break;
			else if(++m > (int) series->l_max)
				m = -(int) series->l_max;
			else if(m == 0)
				break;
			l = abs(m);
		}
	}

	return series;
}
