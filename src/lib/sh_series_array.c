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
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <sphradiometer/math.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                           sh_series_array Object
 *
 * ============================================================================
 */


/*
 * The sh_series_array object is an array of sh_series objects all sharing
 * a single buffer for their coefficients, with each sh_series object's
 * coefficients placed immediately above the coefficients for the preceding
 * sh_series object.  In this way there is a fixed stride from each
 * coefficient to the same coefficient in the next sh_series.  This
 * arrangement allows a sequence of sh_series objects to be operated on
 * collectively using external libraries, for example a sequence of
 * sh_series objects can be Fourier transformed using FFTW by passing the
 * appropriate strides.
 */


/*
 * Create a new sh_series_array object.
 */


struct sh_series_array *sh_series_array_new(int n, unsigned int l_max, int polar)
{
	struct sh_series_array *array = malloc(sizeof(*array));
	struct sh_series *series = malloc(n * sizeof(*series));
	complex double *coeff = malloc(n * sh_series_length(l_max, polar) * sizeof(*coeff));
	int i;

	if(!array || !series || !coeff) {
		free(array);
		free(series);
		free(coeff);
		return NULL;
	}

	for(i = 0; i < n; i++)
		series[i] = (struct sh_series) {
			.l_max = l_max,
			.polar = polar ? 1 : 0,
			.coeff = coeff + i * sh_series_length(l_max, polar)
		};

	array->l_max = l_max;
	array->polar = polar ? 1 : 0;
	array->n = n;
	array->series = series;
	array->coeff = coeff;
	array->stride = sh_series_length(l_max, polar);

	return array;
}


/*
 * Destroy an sh_series_array object.
 */


void sh_series_array_free(struct sh_series_array *array)
{
	if(array) {
		free(array->coeff);
		free(array->series);
	}
	free(array);
}


/*
 * Resize an sh_series_array object.  (only the number of series in the
 * array can be adjusted, not the number of coefficients in each series).
 * If the resize fails, NULL is returned and the original array is left
 * unmodified.
 */


struct sh_series_array *sh_series_array_resize(struct sh_series_array *array, int n)
{
	struct sh_series *series;
	complex double *coeff;
	int i;

	if(n == array->n)
		/* no op */
		return array;

	/* do coeff first because it's bigger, so least likely to work */
	coeff = realloc(array->coeff, n * array->stride * sizeof(*coeff));
	if(!coeff)
		return NULL;
	array->coeff = coeff;

	series = realloc(array->series, n * sizeof(*series));
	if(!series) {
		/* oh oh:  we've successfully reallocated the coefficient
		 * array, but the series reallocation has failed.  Pray we
		 * can unreallocate the coefficients.  If not, oh well. */
		coeff = realloc(array->coeff, array->n * array->stride * sizeof(*coeff));
		if(coeff)
			array->coeff = coeff;
		/* either way, make sure all the series are pointing to the
		 * right place */
		for(i = 0; i < array->n; i++)
			array->series[i].coeff = array->coeff + i * array->stride;
		return NULL;
	}
	array->series = series;

	/* re-initialize all series structures */
	for(i = 0; i < n; i++)
		series[i] = (struct sh_series) {
			.l_max = array->l_max,
			.polar = array->polar,
			.coeff = coeff + i * array->stride,
		};

	array->n = n;

	return array;
}


struct sh_series_array *sh_series_array_resize_zero(struct sh_series_array *array, int n)
{
	const int orig_size = array->n * array->stride;
	const int extra_size = (n - array->n) * array->stride;

	array = sh_series_array_resize(array, n);
	if(!array)
		return NULL;

	if(extra_size > 0)
		memset(array->coeff + orig_size, 0, extra_size * sizeof(*array->coeff));

	return array;
}


/*
 * Set the polar flag of an sh_series_array.  Any new coefficients are zeroed.
 */


struct sh_series_array *sh_series_array_set_polar(struct sh_series_array *array, int polar)
{
	complex double *coeff;
	int i;

	/* sanitize input */
	polar = polar ? 1 : 0;

	/* no-op? */
	if(array->polar == polar)
		return array;

	if(polar) {
		/* making non-azimuthally symmetric sh_series objects
		 * azimuthally symmetric, so the array of coefficients is
		 * getting smaller, so move the coefficients first then
		 * resize the memory */

		/* go forwards, moving the coefficients at the start
		 * downwards to not overwrite anything as we go.  use
		 * sh_series_assign() to do the work of moving the
		 * coefficients */
		for(i = 0; i < array->n; i--) {
			/* create an sh_series object pointing to the new
			 * location of the coefficients */
			struct sh_series series = array->series[i];
			series.polar = polar;
			series.coeff = array->coeff + i * sh_series_length(array->l_max, polar);
			sh_series_assign(&series, &array->series[i]);
			array->series[i].polar = polar;
		}

		/* now resize memory */
		coeff = realloc(array->coeff, array->n * sh_series_length(array->l_max, polar) * sizeof(*coeff));
		if(!coeff)
			return NULL;
		array->coeff = coeff;
	} else {
		/* making azimuthally symmetric sh_series objects
		 * non-azimuthally symmetric, so the array of coefficients is
		 * getting bigger, so resize memory first then move the
		 * coefficients */
		coeff = realloc(array->coeff, array->n * sh_series_length(array->l_max, polar) * sizeof(*coeff));
		if(!coeff)
			return NULL;
		array->coeff = coeff;

		/* go backwards, moving the coefficients at the end upwards
		 * to not overwrite anything as we go.  use
		 * sh_series_assign() to do the work of moving the
		 * coefficients */
		for(i = array->n - 1; i >= 0; i--) {
			/* create an sh_series object pointing to the new
			 * location of the coefficients */
			struct sh_series series = array->series[i];
			series.polar = polar;
			series.coeff = array->coeff + i * sh_series_length(array->l_max, polar);
			sh_series_assign(&series, &array->series[i]);
			array->series[i].polar = polar;
		}
	}

	/* update metadata */
	array->polar = polar;

	/* done */
	return array;
}


/*
 * Make a copy of an sh_series_array object.
 */


struct sh_series_array *sh_series_array_copy(const struct sh_series_array *array)
{
	struct sh_series_array *new = sh_series_array_new(array->n, array->l_max, array->polar);

	if(!new)
		return NULL;

	memcpy(new->coeff, array->coeff, array->n * array->stride * sizeof(*array->coeff));

	return new;
}


/*
 * Assign the coefficients from one sh_series_array object to another.
 * Returns dst on success.  Returns NULL on failure, in which case the
 * coefficients in the destination array are undefined.
 */


struct sh_series_array *sh_series_array_assign(struct sh_series_array *dst, const struct sh_series_array *src)
{
	int i;

	if(dst->n != src->n)
		return NULL;

	for(i = 0; i < dst->n; i++)
		if(!sh_series_assign(&dst->series[i], &src->series[i]))
			return NULL;

	return dst;
}


/*
 * Apply a scale factor to all the elements in an sh_series_array.
 */


struct sh_series_array *sh_series_array_scale(struct sh_series_array *array, complex double z)
{
	int i;

	for(i = 0; i < array->n; i++)
		sh_series_scale(&array->series[i], z);

	return array;
}


/*
 * ============================================================================
 *
 *                    Arithmetic Involving Scalar Vectors
 *
 * ============================================================================
 */


/*
 * Compute the inner product of an array of sh_series objects and a
 * real-valued vector.  The result is an sh_series object equal to
 *
 * 	sum_{i} sh_series[i] * vector[i]
 */


struct sh_series *sh_series_array_dot(struct sh_series *result, const struct sh_series_array *array, const double *vector)
{
	const int n = array->stride;
	const complex double *c = array->coeff;
	const double *last = vector + array->n;

	if((result->l_max != array->l_max) || (result->polar != array->polar))
		return NULL;

	sh_series_zero(result);
	while(vector < last) {
		const double v = *vector++;
		int i;
		for(i = 0; i < n; i++)
			/*result->coeff[i] = cma(*c++, v, result->coeff[i]);*/
			result->coeff[i] += *c++ * v;
	}

	return result;
}


/*
 * Compute the inner product of an array of sh_series objects and a
 * complex-valued vector.  The result is an sh_series object equal to
 *
 * 	sum_{i} sh_series[i] * vector[i]
 */


struct sh_series *sh_series_array_dotc(struct sh_series *result, const struct sh_series_array *array, const complex double *vector)
{
	const int n = array->stride;
	const complex double *c = array->coeff;
	const complex double *last = vector + array->n;

	if((result->l_max != array->l_max) || (result->polar != array->polar))
		return NULL;

	sh_series_zero(result);
	while(vector < last) {
		const complex double v = *vector++;
		int i;
		for(i = 0; i < n; i++)
			/*result->coeff[i] = ccma(*c++, v, result->coeff[i]);*/
			result->coeff[i] += *c++ * v;
	}

	return result;
}


/*
 * Apply a real-valued window function vector to the sh_series objects in
 * an sh_series_array object.
 */


struct sh_series_array *sh_series_array_window(struct sh_series_array *array, const double *window)
{
	const int n = array->stride;
	complex double *c = array->coeff;
	const double *last_w = window + array->n;

	while(window < last_w) {
		const complex double *last_c = c + n;
		const double w = *window++;
		while(c < last_c)
			*c++ *= w;
	}

	return array;
}


/*
 * Apply a complex-valued window function vector to the sh_series objects
 * in an sh_series_array object.
 */


struct sh_series_array *sh_series_array_windowc(struct sh_series_array *array, const complex double *window)
{
	const int n = array->stride;
	complex double *c = array->coeff;
	const complex double *last_w = window + array->n;

	while(window < last_w) {
		const complex double *last_c = c + n;
		const complex double w = *window++;
		while(c < last_c)
			*c++ *= w;
	}

	return array;
}


/*
 * ============================================================================
 *
 *                             Fourier Transforms
 *
 * ============================================================================
 */


/*
 * Compute the forward Fourier transform of an array of sh_series objects.
 * The vector of c_{0,0} coefficients is replaced by its Fourier transform,
 * the vector c_{1,-1} coefficients is replaced with its Fourier transform,
 * and so on.  Following the transform, sh_series object #0 in the
 * sh_series_array contains the DC components for all coefficients, and the
 * rest are stored in order in the other sh_series objects.
 *
 * NOTE:  it is assumed this, and the reverse transforms, are very
 * infrequent procedures.  If not, new functions should be introduced to
 * precompute the plans and re-use them.
 */


struct sh_series_array *sh_series_array_forward_fft(struct sh_series_array *array)
{
	const int n[] = {array->n};
	/* have to use FFTW_ESTIMATE or the array's contents are destroyed */
	fftw_plan plan = fftw_plan_many_dft(1, n, array->stride, array->coeff, NULL, array->stride, 1, array->coeff, NULL, array->stride, 1, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(plan);

	fftw_destroy_plan(plan);

	return array;
}


/*
 * Compute the inverse Fourier transform of an array of sh_series objects.
 * On input, sh_series object #0 in the sh_series_array contains the DC
 * components.
 */


struct sh_series_array *sh_series_array_reverse_fft(struct sh_series_array *array)
{
	const int n[] = {array->n};
	/* have to use FFTW_ESTIMATE or the array's contents are destroyed */
	fftw_plan plan = fftw_plan_many_dft(1, n, array->stride, array->coeff, NULL, array->stride, 1, array->coeff, NULL, array->stride, 1, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftw_execute(plan);

	fftw_destroy_plan(plan);
	return array;
}
