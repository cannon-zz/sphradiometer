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
#include <stdlib.h>
#include <string.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                              sh_series Object
 *
 * ============================================================================
 */


/*
 * Create a new sh_series object, with maximum l of l_max.
 */


struct sh_series *sh_series_new(unsigned int l_max, int polar)
{
	struct sh_series *series = malloc(sizeof(*series));
	complex double *coeff = malloc(sh_series_length(l_max, polar) * sizeof(*coeff));

	if(!series || !coeff) {
		free(series);
		free(coeff);
		return NULL;
	}

	series->l_max = l_max;
	series->polar = polar ? 1 : 0;
	series->coeff = coeff;

	return series;
}


/*
 * Create a new sh_series object, with maximum l of l_max, and the
 * coefficients zeroed.  This is marginally faster than sh_series_new() +
 * sh_series_zero().
 */


struct sh_series *sh_series_new_zero(unsigned int l_max, int polar)
{
	struct sh_series *series = malloc(sizeof(*series));
	complex double *coeff = calloc(sh_series_length(l_max, polar), sizeof(*coeff));

	if(!series || !coeff) {
		free(series);
		free(coeff);
		return NULL;
	}

	series->l_max = l_max;
	series->polar = polar ? 1 : 0;
	series->coeff = coeff;

	return series;
}


/*
 * Create a copy of an sh_series object.
 */


struct sh_series *sh_series_copy(const struct sh_series *series)
{
	struct sh_series *new = sh_series_new(series->l_max, series->polar);

	if(!new)
		return NULL;

	memcpy(new->coeff, series->coeff, sh_series_length(series->l_max, series->polar) * sizeof(*new->coeff));

	return new;
}


/*
 * Assign the coefficients from one sh_series object to another.  Returns
 * dst on success, NULL on failure.  This operation is permitted to be
 * lossy:  the destination is allowed to possess a proper subset of the
 * source's coefficients.  If calling code requires error checking for
 * this, it must provide it itself.
 */


struct sh_series *sh_series_assign(struct sh_series *dst, const struct sh_series *src)
{
	/* max m that needs to be set in the target */
	int dst_m_max = dst->polar ? 0 : dst->l_max;
	/* max m to copy from src, the rest are set to 0 */
	int m_max = (src->polar || dst->polar) ? 0 : (src->l_max < dst->l_max) ? src->l_max : dst->l_max;
	int m;

	/* NOTE:  this algorithm is a little more complicated than it needs
	 * to be because it works even in some situations in which the
	 * coefficient arrays of dst and src overlap.  this is done to
	 * allow this function to be used internally to help implement some
	 * of the other operations of the library, but it does not work in
	 * all such situations and is not documented behaviour and must not
	 * be relied on by external code */

	if(dst->l_max < src->l_max) {
		/* move through the coefficient array in the forward
		 * direction, copying one set of m=constant coefficients at
		 * a time.  zero any extra coefficients */
		for(m = 0; m <= m_max; m++)
			memmove(dst->coeff + sh_series_moffset(dst->l_max, m), src->coeff + sh_series_moffset(src->l_max, m), (dst->l_max + 1 - m) * sizeof(*dst->coeff));
		for(; m <= dst_m_max; m++)
			memset(dst->coeff + sh_series_moffset(dst->l_max, m), 0, (dst->l_max + 1 - m) * sizeof(*dst->coeff));
		for(m = -dst_m_max; m < -m_max; m++)
			memset(dst->coeff + sh_series_moffset(dst->l_max, m), 0, (dst->l_max + 1 + m) * sizeof(*dst->coeff));
		for(m = -m_max; m < 0; m++)
			memmove(dst->coeff + sh_series_moffset(dst->l_max, m), src->coeff + sh_series_moffset(src->l_max, m), (dst->l_max + 1 + m) * sizeof(*dst->coeff));
	} else if(dst->l_max == src->l_max) {
		if(dst->polar || !src->polar) {
			/* either src or dst have identical parameters, or
			 * dst has azimuthal symmetry and src doesn't,
			 * either way dst's coefficients are a subset of
			 * src's so we can simply copy them verbatim */
			memmove(dst->coeff, src->coeff, sh_series_length(dst->l_max, dst->polar) * sizeof(*dst->coeff));
		} else {
			/* src has azimuthal symmetry but dst doesn't, so
			 * copy the coefficients that are available in src,
			 * then zero the rest of dst's coefficients. */
			memmove(dst->coeff, src->coeff, sh_series_length(src->l_max, src->polar) * sizeof(*dst->coeff));
			memset(dst->coeff + sh_series_length(src->l_max, src->polar), 0, (sh_series_length(dst->l_max, dst->polar) - sh_series_length(src->l_max, src->polar)) * sizeof(*dst->coeff));
		}
	} else /* (dst->l_max > src->l_max) */ {
		/* move through the coefficient array in the reverse
		 * direction, copying one set of m=constant coefficients at
		 * a time */
		for(m = -1; m >= -m_max; m--) {
			memmove(dst->coeff + sh_series_moffset(dst->l_max, m), src->coeff + sh_series_moffset(src->l_max, m), (src->l_max + 1 + m) * sizeof(*dst->coeff));
			memset(dst->coeff + sh_series_moffset(dst->l_max, m) + (src->l_max + 1 + m), 0, (dst->l_max - src->l_max) * sizeof(*dst->coeff));
		}
		for(; m >= -dst_m_max; m--)
			memset(dst->coeff + sh_series_moffset(dst->l_max, m), 0, (dst->l_max + 1 + m) * sizeof(*dst->coeff));
		for(m = dst_m_max; m > m_max; m--)
			memset(dst->coeff + sh_series_moffset(dst->l_max, m), 0, (dst->l_max + 1 - m) * sizeof(*dst->coeff));
		for(m = m_max; m >= 0; m--) {
			memmove(dst->coeff + sh_series_moffset(dst->l_max, m), src->coeff + sh_series_moffset(src->l_max, m), (src->l_max + 1 - m) * sizeof(*dst->coeff));
			memset(dst->coeff + sh_series_moffset(dst->l_max, m) + (src->l_max + 1 - m), 0, (dst->l_max - src->l_max) * sizeof(*dst->coeff));
		}
	}

	return dst;
}


/*
 * Resize an sh_series object to the new maximum l, preserving the function
 * represented by the coefficients (up to the smaller of the new and old
 * harmonic orders).  Returns NULL if the resize failed.  The resize is
 * done in place.  If the resize fails, then if the series was being shrunk
 * in size then the contents are undefined, otherwise they are unchanged.
 */


struct sh_series *sh_series_resize(struct sh_series *series, unsigned int l_max)
{
	/* clone series to use as a source, adjust the harmonic order of
	 * series to use as a target, then use sh_series_assign() from the
	 * source to the target to shuffle and initialize coefficients as
	 * needed. */

	struct sh_series wrapper = *series;
	complex double *coeff;

	series->l_max = l_max;

	if(wrapper.l_max > l_max) {
		/* array is getting smaller:  realloc after move */
		if(!sh_series_assign(series, &wrapper))
			return NULL;
		coeff = realloc(series->coeff, sh_series_length(series->l_max, series->polar) * sizeof(*coeff));
		if(!coeff) {
			/* ignore failure */
		} else
			series->coeff = coeff;
		return series;
	} else if(wrapper.l_max == l_max) {
		/* no-op? */
		return series;
	} else {	/* wrapper.l_max < l_max */
		/* array is getting bigger:  realloc before move */
		coeff = realloc(series->coeff, sh_series_length(series->l_max, series->polar) * sizeof(*coeff));
		if(!coeff)
			return NULL;
		wrapper.coeff = series->coeff = coeff;
		return sh_series_assign(series, &wrapper);
	}
}


/*
 * Zero the coefficients of an sh_series object.
 */


struct sh_series *sh_series_zero(struct sh_series *series)
{
	memset(series->coeff, 0, sh_series_length(series->l_max, series->polar) * sizeof(*series->coeff));

	return series;
}


/*
 * Set the polar flag of an sh_series.  Any new coefficients are zeroed.
 */


struct sh_series *sh_series_set_polar(struct sh_series *series, int polar)
{
	complex double *coeff;

	/* sanitize input */
	polar = polar ? 1 : 0;

	/* no-op? */
	if(series->polar == polar)
		return series;

	/* resize the coefficients array */
	coeff = realloc(series->coeff, sh_series_length(series->l_max, polar) * sizeof(*coeff));
	if(!coeff)
		return NULL;
	series->coeff = coeff;

	if(!polar)
		/* array got bigger.  zero the new coefficients */
		memset(series->coeff + sh_series_length(series->l_max, series->polar), 0, (sh_series_length(series->l_max, polar) - sh_series_length(series->l_max, series->polar)) * sizeof(*series->coeff));

	/* update metadata.  done */
	series->polar = polar;
	return series;
}


/*
 * Free an sh_series object.
 */


void sh_series_free(struct sh_series *series)
{
	if(series)
		free(series->coeff);
	free(series);
}


/*
 * Return the coefficient for the given l and m.
 */


complex double sh_series_get(const struct sh_series *series, unsigned int l, int m)
{
	return series->coeff[sh_series_lmoffset(series, l, m)];
}


/*
 * Set the coefficient for the given l and m.
 */


complex double sh_series_set(struct sh_series *series, unsigned int l, int m, complex double val)
{
	return series->coeff[sh_series_lmoffset(series, l, m)] = val;
}


/*
 * Are any coefficients NaN?
 */


int sh_series_is_nan(const struct sh_series *series)
{
	complex double *coeff = series->coeff;
	complex double *last = coeff + sh_series_length(series->l_max, series->polar);
	for(; coeff < last; coeff++)
		if(isnan(creal(*coeff)) || isnan(cimag(*coeff)))
			return 1;
	return 0;
}


/*
 * Are any coefficients inf?
 */


int sh_series_is_inf(const struct sh_series *series)
{
	complex double *coeff = series->coeff;
	complex double *last = coeff + sh_series_length(series->l_max, series->polar);
	for(; coeff < last; coeff++)
		if(isinf(creal(*coeff)) || isinf(cimag(*coeff)))
			return 1;
	return 0;
}


/*
 * ============================================================================
 *
 *                     Simple sh_series Object Arithmetic
 *
 * ============================================================================
 */


/*
 * Add z * sh_series object b to sh_series object a in place, returning a
 * pointer to a or NULL on failure.  "a" and "b" need not have the same
 * order nor the same azimuthal symmetry, but the operation is not
 * permitted to be lossy so "a" cannot have a subset of the coefficients of
 * "b".
 */


struct sh_series *sh_series_add(struct sh_series *a, const complex double z, const struct sh_series *b)
{
	const int m_max = b->polar ? 0 : b->l_max;
	int m;

	if((a->l_max == b->l_max) && (a->polar == b->polar)) {
		complex double *dst = a->coeff;
		complex double *src = b->coeff;
		complex double *last = src + sh_series_length(b->l_max, b->polar);
		while(src < last)
			*dst++ += z * *src++;
		return a;
	}

	if((a->l_max < b->l_max) || (a->polar && !b->polar))
		return NULL;

	for(m = -m_max; m <= m_max; m++) {
		complex double *dst = a->coeff + sh_series_moffset(a->l_max, m);
		complex double *src = b->coeff + sh_series_moffset(b->l_max, m);
		complex double *last = src + b->l_max + 1 - abs(m);
		while(src < last)
			*dst++ += z * *src++;
	}

	return a;
}


/*
 * Multiply the coefficients in the sh_series object a by the complex
 * factor z, returning a pointer to a or NULL on failure.
 */


struct sh_series *sh_series_scale(struct sh_series *a, complex double z)
{
	complex double *c = a->coeff;
	const complex double *c_last = c + sh_series_length(a->l_max, a->polar);

	while(c < c_last)
		*c++ *= z;

	return a;
}


/*
 * Replace an sh_series object with its complex conjugate.
 */


struct sh_series *sh_series_conj(struct sh_series *a)
{
	complex double *c = a->coeff;
	const complex double *c_last = c + sh_series_length(a->l_max, a->polar);

	for(; c < c_last; c++)
		*c = conj(*c);

	return a;
}


/*
 * Set to zero any coefficients whose magnitude as a fraction of the
 * largest magnitude coefficient is less than epsilon.
 */


struct sh_series *sh_series_clip(struct sh_series *series, double epsilon)
{
	complex double *c;
	const complex double *c_last = series->coeff + sh_series_length(series->l_max, series->polar);
	double threshold = 0.0;

	/* find the magnitude of the largest coefficient */
	for(c = series->coeff; c < c_last;)
		threshold = fmax(threshold, cabs(*c++));

	/* no-op? */
	if(threshold == 0.0)
		return series;

	/* convert fractional threshold to absolute */
	threshold *= fabs(epsilon);

	/* clip values */
	for(c = series->coeff; c < c_last; c++)
		if(cabs(*c) < threshold)
			*c = 0.0;

	return series;
}


/*
 * Compute the inner product of two functions on the sphere.  This uses an
 * algorithm that retains high accuracy when summing large numbers of
 * coefficients, but is slower than a simple sum-in-a-loop, and requires
 * memory for scratch space.  Returns the inner product on success, and
 * complex NaN on failure.
 */


static int _complex_double_cmp(const void *a, const void *b)
{
	const double z = cabs(*(complex double *) a - *(complex double *) b);

	return z > 0.0 ? +1 : z < 0.0 ? -1 : 0;
}


complex double sh_series_dot(const struct sh_series *a, const struct sh_series *b)
{
	const int n = sh_series_length(a->l_max, a->polar);
	complex double *tmp = malloc(n * sizeof(*tmp));
	const complex double *a_coeff = a->coeff + n;
	const complex double *b_coeff = b->coeff + n;
	complex double c;
	int i;

	if(!tmp || (b->l_max != a->l_max) || (b->polar != a->polar)) {
		/* FIXME: add support for these cases */
		free(tmp);
		return nan("") + I * nan("");
	}

	/* put the products of the highest-order coefficients into the
	 * lowest addresses in tmp[] in the belief that this makes tmp
	 * mostly sorted in smallest to largest order */
	for(i = 0; i < n; i++)
		tmp[i] = *--a_coeff * conj(*--b_coeff);

	/* sort from smallest magnitude to largest magnitude */
	qsort(tmp, n, sizeof(*tmp), _complex_double_cmp);

	/* sum */
	c = 0.0;
	for(i = 0; i < n; i++)
		c += tmp[i];

	free(tmp);
	return c;
}


/*
 * ============================================================================
 *
 *                              Differentiation
 *
 * ============================================================================
 */


/*
 * Replace an sh_series object with its Laplacian.  Spherical harmonics are
 * eigenfunctions of the Laplacian,
 *
 * 	\grad^{2} Y_{lm}(\hat{s}) = -l (l + 1) Y_{lm}(\hat{s})
 *
 * so all we do is run down the list of coefficients multiplying each by
 * the eigenvalue.
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
