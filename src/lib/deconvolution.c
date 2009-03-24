/*
 * Copyright (C) 2006  Kipp C. Cannon
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
#include <stdlib.h>
#include <radiometer/sh_series.h>
#include <radiometer/deconvolution.h>


/*
 * ============================================================================
 *
 *                            Image Deconvolution
 *
 * ============================================================================
 */


/*
 * Starting from an sh_series_array object containing n
 * linearly-independant functions on the sphere, use the Gram-Schmidt
 * process to construct a set of orthonormal basis functions complete to
 * order l_max of the array.  The first n basis functions span the space
 * spanned by the n original functions.  The basis functions are returned
 * in the sh_series_array, whose size may be increased to accomodate them,
 * and the function's return value is the vector of coefficients defining
 * the sum of the original n functions in terms of the new orthonormal
 * basis functions.
 */


complex double *sh_series_array_orthogonalize(struct sh_series_array *array)
{
	const int n = array->n;
	complex double *vector = calloc(array->stride, sizeof(*vector));
	int i, j;

	if(!vector || (array->n > array->stride)) {
		free(vector);
		return NULL;
	}

	/* make the sh_series_array object square (give it as many elements
	 * as each has coefficients) */
	if(!sh_series_array_resize_zero(array, array->stride)) {
		free(vector);
		return NULL;
	}

	/* fill unused parts of array with new functions that are linearly
	 * independant of the first n.  FIXME:  although this hack almost
	 * certainly works (what are the odds of one of the baselines
	 * producing a brightness distribution exactly equal to one of the
	 * Y_{lm}'s?), it would be nice to be certain. */
	for(i = n; i < array->n; i++)
		array->series[i].coeff[i] = 1.0;

	/* Gram-Schmidt process */
	for(i = 0; i < array->n; i++) {
		complex double c;
		for(j = 0; j < i; j++) {
			c = sh_series_dot(&array->series[j], &array->series[i]);
			if(i < n)
				vector[j] += c;
			sh_series_sub(&array->series[i], c, &array->series[j]);
		}
		c = csqrt(sh_series_dot(&array->series[i], &array->series[i]));
		if(i < n)
			vector[i] += c;
		sh_series_scale(&array->series[i], 1.0 / c);
	}

	return vector;
}
