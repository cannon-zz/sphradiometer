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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                       Spherical Harmonic Evaluation
 *
 * ============================================================================
 */


/*
 * Return Y_{lm}(theta, phi), and Y^{*}_{lm}(theta, phi).
 */


complex double sh_series_Y(unsigned int l, int m, double theta, double phi)
{
	return ((m < 0 && m & 1) ? -1.0 : +1.0) * (cos(m * phi) + I * sin(m * phi)) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
}


complex double sh_series_Yconj(unsigned int l, int m, double theta, double phi)
{
	return ((m < 0 && m & 1) ? -1.0 : +1.0) * (cos(m * phi) - I * sin(m * phi)) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
}


/*
 * Compute an array of Y_{lm}(theta, phi), and Y^{*}_{lm}(theta, phi) for l
 * in [|m|, ..., l_max].
 */


complex double *sh_series_Y_array(complex double *array, unsigned int l_max, int m, double theta, double phi)
{
	double tmp[l_max + 1 - abs(m)];
	const complex double factor = ((m < 0) && (m & 1) ? -1.0 : +1.0) * cexp(I * m * phi);
	int i;

	gsl_sf_legendre_sphPlm_array(l_max, abs(m), cos(theta), tmp);

	for(i = 0; i <= (int) l_max - abs(m); i++)
		array[i] = factor * tmp[i];

	return array;
}


complex double *sh_series_Yconj_array(complex double *array, unsigned int l_max, int m, double theta, double phi)
{
	double tmp[l_max + 1 - abs(m)];
	const complex double factor = ((m < 0) && (m & 1) ? -1.0 : +1.0) * cexp(I * -m * phi);
	int i;

	gsl_sf_legendre_sphPlm_array(l_max, abs(m), cos(theta), tmp);

	for(i = 0; i <= (int) l_max - abs(m); i++)
		array[i] = factor * tmp[i];

	return array;
}


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
 * coefficients zeroed.  This marginally faster than sh_series_new() +
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
 * dst on success, NULL on failure.
 */


struct sh_series *sh_series_assign(struct sh_series *dst, const struct sh_series *src)
{
	if(dst->l_max != src->l_max)
		return NULL;

	if(dst->polar || !src->polar) {
		/* either src or dst have identical parameters, or dst has
		 * azimuthal symmetry and src doesn't, either way dst's
		 * coefficients are a subset of src's so we can simply copy
		 * them verbatim */
		memcpy(dst->coeff, src->coeff, sh_series_length(dst->l_max, dst->polar) * sizeof(*dst->coeff));
	} else {
		/* src has azimuthal symmetry but dst doesn't, so zero all
		 * of dst's coefficients, then copy the ones that are set
		 * in src */
		memset(dst->coeff, 0, sh_series_length(dst->l_max, dst->polar) * sizeof(*dst->coeff));
		memcpy(dst->coeff, src->coeff, sh_series_length(src->l_max, src->polar) * sizeof(*dst->coeff));
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


static struct sh_series *_sh_series_resize(struct sh_series *series, unsigned int l_max)
{
	complex double *newcoeff;
	int m;
	/* to what order the coefficients must be preserved (smaller of the
	 * new and old orders) */
	const int preserve_l_max = series->l_max < l_max ? series->l_max : l_max;
	/* corresponding m limit */
	const int m_max = series->polar ? 0 : preserve_l_max;

	if(l_max < series->l_max) {
		/* shrinking:  move coefficients before resize */

		/* copy the coefficients, moving foward through the array
		 * */
		for(m = 0; m <= m_max; m++)
			memcpy(series->coeff + sh_series_moffset(l_max, m), series->coeff + sh_series_moffset(series->l_max, m), (preserve_l_max - abs(m) + 1) * sizeof(*series->coeff));
		for(m = -m_max; m < 0; m++)
			memcpy(series->coeff + sh_series_moffset(l_max, m), series->coeff + sh_series_moffset(series->l_max, m), (preserve_l_max - abs(m) + 1) * sizeof(*series->coeff));

		/* reallocate coefficient array */
		newcoeff = realloc(series->coeff, sh_series_length(l_max, series->polar) * sizeof(*series->coeff));
		if(!newcoeff)
			return NULL;
		series->coeff = newcoeff;
	} else {
		/* not shrinking:  move coefficients after resize */

		/* reallocate coefficient array */
		newcoeff = realloc(series->coeff, sh_series_length(l_max, series->polar) * sizeof(*series->coeff));
		if(!newcoeff)
			return NULL;
		series->coeff = newcoeff;

		/* copy the coefficients, moving backward through the
		 * array, and zeroing new coefficients */
		for(m = -1; m >= -m_max; m--) {
			memset(series->coeff + sh_series_params_lmoffset(series->l_max, preserve_l_max + 1, m), 0, (l_max - preserve_l_max) * sizeof(*series->coeff));
			memcpy(series->coeff + sh_series_moffset(l_max, m), series->coeff + sh_series_moffset(series->l_max, m), (preserve_l_max - abs(m) + 1) * sizeof(*series->coeff));
		}
		for(m = m_max; m >= 0; m--) {
			memset(series->coeff + sh_series_params_lmoffset(series->l_max, preserve_l_max + 1, m), 0, (l_max - preserve_l_max) * sizeof(*series->coeff));
			memcpy(series->coeff + sh_series_moffset(l_max, m), series->coeff + sh_series_moffset(series->l_max, m), (preserve_l_max - abs(m) + 1) * sizeof(*series->coeff));
		}
	}

	series->l_max = l_max;

	return series;
}


struct sh_series *sh_series_resize(struct sh_series *series, unsigned int l_max)
{
	/* no-op? */
	if(series->l_max == l_max)
		return series;

	/* do the resize */
	return _sh_series_resize(series, l_max);
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
	/* sanitize input */
	polar = polar ? 1 : 0;

	/* no-op? */
	if(series->polar == polar)
		return series;

	/* use _sh_series_resize() do to the work */
	series->polar = polar;
	return _sh_series_resize(series, series->l_max);
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
 * ============================================================================
 *
 *                     Simple sh_series Object Arithmetic
 *
 * ============================================================================
 */


/*
 * Add z * sh_series object b to sh_series object a in place, returning a
 * pointer to a or NULL on failure.
 */


struct sh_series *sh_series_add(struct sh_series *a, const complex double z, const struct sh_series *b)
{
	complex double *dst = a->coeff;
	complex double *src = b->coeff;
	unsigned int n = sh_series_length(b->l_max, b->polar);

	if((a->l_max != b->l_max) || (a->polar != b->polar))
		return NULL;

	while(n--)
		*dst++ += z * *src++;

	return a;
}


/*
 * Add z * sh_series object b to sh_series object a in place, returning a
 * pointer to a or NULL on failure.  Unlike the _add() function which
 * requires "a" and "b" to have identical orders, this function allows "a"
 * to have a higher order than "b".  For being more general, this
 * function is slower than _add().
 */


struct sh_series *sh_series_add_into(struct sh_series *a, const complex double z, const struct sh_series *b)
{
	const int m_max = b->polar ? 0 : b->l_max;
	int l, m;

	if((a->l_max < b->l_max) || (a->polar && !b->polar))
		return NULL;

	for(m = -m_max; m <= m_max; m++) {
		complex double *dst = a->coeff + sh_series_moffset(a->l_max, m);
		complex double *src = b->coeff + sh_series_moffset(b->l_max, m);
		for(l = abs(m); l <= (int) b->l_max; l++)
			*dst++ += z * *src++;
	}

	return a;
}


/*
 * Subtract z *  sh_series object b from sh_series object a in place,
 * returning a pointer to a or NULL on failure.
 */


struct sh_series *sh_series_sub(struct sh_series *a, const complex double z, const struct sh_series *b)
{
	complex double *dst = a->coeff;
	complex double *src = b->coeff;
	unsigned int n = sh_series_length(b->l_max, b->polar);

	if((a->l_max != b->l_max) || (a->polar != b->polar))
		return NULL;

	while(n--)
		*dst++ -= z * *src++;

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
 *                       Spherical Harmonic Transforms
 *
 * ============================================================================
 */


/*
 * Evaluate an sh_series object at theta, phi.
 */


complex double sh_series_eval(const struct sh_series *series, double theta, double phi)
{
	complex double vals[series->l_max + 1];
	complex double *coeff = series->coeff;
	int l, m;
	complex double val = 0.0;

	sh_series_Y_array(vals, series->l_max, 0, theta, phi);
	for(l = 0; l <= (int) series->l_max; l++)
		val += *(coeff++) * vals[l];
	if(!series->polar) {
		for(m = 1; m <= (int) series->l_max; m++) {
			complex double *v = vals;
			sh_series_Y_array(v, series->l_max, m, theta, phi);
			for(l = m; l <= (int) series->l_max; l++)
				val += *(coeff++) * *(v++);
		}
		for(m = -(int) series->l_max; m < 0; m++) {
			complex double *v = vals;
			sh_series_Y_array(v, series->l_max, m, theta, phi);
			for(l = -m; l <= (int) series->l_max; l++)
				val += *(coeff++) * *(v++);
		}
	}

	return val;
}


/*
 * Given a function of theta, phi, evaluate it on a 2-dimensional mesh of
 * points, for later input to a projection function.
 */


static void samples_from_l_max(unsigned l_max, int *ntheta, int *nphi)
{
	*ntheta = 32 * l_max;
	*nphi = 2 * l_max;
}


complex double *sh_series_mesh_new(unsigned int l_max, int *ntheta, int *nphi)
{
	int nt, np;
	complex double *f;

	samples_from_l_max(l_max, &nt, &np);

	f = malloc(nt * np * sizeof(*f));

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	return f;
}


double *sh_series_real_mesh_new(unsigned int l_max, int *ntheta, int *nphi)
{
	int nt, np;
	double *f;

	samples_from_l_max(l_max, &nt, &np);

	f = malloc(nt * np * sizeof(*f));

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	return f;
}


/*
 * This fills the array with samples of f in the same order and evaluated
 * at the same points on the sphere as required by the functions in S2kit.
 */


complex double *sh_series_mesh_from_func(unsigned int l_max, complex double (*func)(double, double, void *), void *data, int *ntheta, int *nphi)
{
	int nt, np;
	complex double *f = sh_series_mesh_new(l_max, &nt, &np);
	const double dtheta = M_PI / nt;
	const double dphi = 2 * M_PI / np;
	int i, j;

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	if(f)
		/* fill f[][] with samples from the function to be projected */
		for(i = 0; i < nt; i++) {
			const double theta = dtheta * (i + 0.5);
			for(j = 0; j < np; j++)
				*(f + i * np + j) = func(theta, dphi * j, data);
		}

	return f;
}


double *sh_series_mesh_from_realfunc(unsigned int l_max, double (*func)(double, double, void *), void *data, int *ntheta, int *nphi)
{
	int nt, np;
	double *f = sh_series_real_mesh_new(l_max, &nt, &np);
	const double dtheta = M_PI / nt;
	const double dphi = 2 * M_PI / np;
	int i, j;

	if(ntheta)
		*ntheta = nt;
	if(nphi)
		*nphi = np;

	if(f)
		/* fill f[][] with samples from the function to be projected */
		for(i = 0; i < nt; i++) {
			const double theta = dtheta * (i + 0.5);
			for(j = 0; j < np; j++)
				*(f + i * np + j) = func(theta, dphi * j, data);
		}

	return f;
}


/*
 * Given a 2-D mesh of samples, project the data onto spherical harmonics,
 * and return the coefficients as an sh_series object.  Returns NULL on
 * failure.  See the functions above for how the mesh of samples is
 * defined.  Note that, as a side effect, these functions modify the
 * contents of the mesh array.
 */


struct sh_series *sh_series_from_mesh(struct sh_series *series, complex double *mesh)
{
	int ntheta, nphi;
	double dtheta, dphi;
	int n[1];
	complex double *F;
	double *P;
	fftw_plan plan;
	int l, m, i, j;

	samples_from_l_max(series->l_max, &ntheta, &nphi);
	dtheta = M_PI / ntheta;
	dphi = 2 * M_PI / nphi;
	n[0] = nphi;
	F = malloc(ntheta * nphi * sizeof(*F));
	P = malloc(ntheta * (series->l_max + 1) * sizeof(*P));

	if(!F || !P) {
		free(F);
		free(P);
		return NULL;
	}

	/* construct the FFTW plan */
	plan = fftw_plan_many_dft(1, n, ntheta, mesh, NULL, 1, nphi, F, NULL, 1, nphi, FFTW_FORWARD, FFTW_ESTIMATE);

	/* multiply the function samples by dOmega = sin theta dtheta dphi
	 * to save having to multiply by these factors in the inner-most
	 * loop below */
	for(i = 0; i < ntheta; i++) {
		const double domega = sin(dtheta * (i + 0.5)) * dtheta * dphi;
		for(j = 0; j < nphi; j++)
			*(mesh + i * nphi + j) *= domega;
	}

	/* execute the FFTW plan, after this F[][] contains the
	 * H_{m}(theta) sin(theta) dtheta dphi factors */
	fftw_execute(plan);

	/* for each non-negative m, */
	for(m = 0; m <= (int) (series->polar ? 0 : series->l_max); m++) {
		/* populate P[][] with (normalized) P_{l m}(cos theta) */
		for(i = 0; i < ntheta; i++)
			gsl_sf_legendre_sphPlm_array(series->l_max, m, cos(dtheta * (i + 0.5)), P + i * (series->l_max + 1) + m);
		/* for each l s.t. m <= l <= l_max, */
		for(l = m; l <= (int) series->l_max; l++) {
			/* integrate H_{m}(theta) sin(theta) P_{l m}(cos
			 * theta) over theta (+m case) */
			complex double c = 0.0;
			for(i = 0; i < ntheta; i++)
				c += *(F + i * nphi + m) * *(P + i * (series->l_max + 1) + l);
			sh_series_set(series, l, m, c);
		}
		if(m) {
			const double sign = m & 1 ? -1.0 : +1.0;
			/* for each l s.t. m <= l <= l_max, */
			for(l = m; l <= (int) series->l_max; l++) {
				/* integrate H_{-m}(theta) sin(theta) P_{l
				 * -m}(cos theta) over theta (-m case) */
				complex double c = 0.0;
				for(i = 0; i < ntheta; i++)
					c += *(F + i * nphi + (nphi - m)) * *(P + i * (series->l_max + 1) + l);
				sh_series_set(series, l, -m, sign * c);
			}
		}
	}

	/* clean up */
	fftw_destroy_plan(plan);
	free(P);
	free(F);

	return sh_series_clip(series, 1e-10);
}


struct sh_series *sh_series_from_realmesh(struct sh_series *series, double *mesh)
{
	int ntheta, nphi;
	double dtheta, dphi;
	int n[1];
	complex double *F;
	double *P;
	fftw_plan plan;
	int l, m, i, j;

	samples_from_l_max(series->l_max, &ntheta, &nphi);
	dtheta = M_PI / ntheta;
	dphi = 2 * M_PI / nphi;
	n[0] = nphi;
	F = malloc(ntheta * (nphi / 2 + 1) * sizeof(*F));
	P = malloc(ntheta * (series->l_max + 1) * sizeof(*P));

	if(!F || !P) {
		free(F);
		free(P);
		return NULL;
	}

	/* construct the FFTW plan */
	plan = fftw_plan_many_dft_r2c(1, n, ntheta, mesh, NULL, 1, nphi, F, NULL, 1, nphi / 2 + 1, FFTW_ESTIMATE);

	/* multiply the function samples by dOmega = sin theta dtheta dphi
	 * to save having to multiply by these factors in the inner-most
	 * loop below */
	for(i = 0; i < ntheta; i++) {
		const double domega = sin(dtheta * (i + 0.5)) * dtheta * dphi;
		for(j = 0; j < nphi; j++)
			*(mesh + i * nphi + j) *= domega;
	}

	/* execute the FFTW plan, after this F[][] contains the
	 * H_{m}(theta) sin(theta) dtheta dphi factors */
	fftw_execute(plan);

	/* for each non-negative m, */
	for(m = 0; m <= (int) (series->polar ? 0 : series->l_max); m++) {
		/* populate P[][] with (normalized) P_{lm}(cos theta) */
		for(i = 0; i < ntheta; i++)
			gsl_sf_legendre_sphPlm_array(series->l_max, m, cos(dtheta * (i + 0.5)), P + i * (series->l_max + 1) + m);
		/* for each l s.t. m <= l <= l_max, */
		for(l = m; l <= (int) series->l_max; l++) {
			/* integrate H_{m}(theta) sin(theta) P_{lm}(cos
			 * theta) over theta for +/-m */
			complex double c = 0.0;
			for(i = 0; i < ntheta; i++)
				c += *(F + i * (series->l_max + 1) + m) * *(P + i * (series->l_max + 1) + l);
			sh_series_set(series, l, -m, (m & 1 ? -1 : 1) * conj(sh_series_set(series, l, m, c)));
		}
	}

	/* clean up */
	fftw_destroy_plan(plan);
	free(P);
	free(F);

	return sh_series_clip(series, 1e-10);
}


/*
 * Given a function of theta, phi, project the function onto spherical
 * harmonics upto and including l = l_max, and return the coefficients as
 * an sh_series object.  Returns NULL on failure.  (Convenience wrappers of
 * sh_series_mesh_from_func() and sh_series_from_mesh()).
 */


struct sh_series *sh_series_from_func(struct sh_series *series, complex double (*func)(double, double, void *), void *data)
{
	complex double *mesh = sh_series_mesh_from_func(series->l_max, func, data, NULL, NULL);

	if(!mesh)
		return NULL;

	series = sh_series_from_mesh(series, mesh);
	free(mesh);

	return series;
}


struct sh_series *sh_series_from_realfunc(struct sh_series *series, double (*func)(double, double, void *), void *data)
{
	double *mesh = sh_series_mesh_from_realfunc(series->l_max, func, data, NULL, NULL);

	if(!mesh)
		return NULL;

	series = sh_series_from_realmesh(series, mesh);
	free(mesh);
	return series;
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


/*
 * ============================================================================
 *
 *                                  File I/O
 *
 * ============================================================================
 */


/*
 * Print the non-zero coefficients in an sh_series object.
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


/*
 * Write an sh_series object to a file.  Warning:  the data is not
 * portable;  an sh_series object can only be reconstructed from the file
 * on a machine with the same word size, alignment, and endianness.
 * Returns 0 on success, -1 on failure (use ferror() or errno for details).
 *
 * This is a hack, don't use in production code.
 */


int sh_series_write(const struct sh_series *series, FILE *file)
{
	size_t n;

	n = fwrite(&series->l_max, sizeof(series->l_max), 1, file);
	n += fwrite(&series->polar, sizeof(series->polar), 1, file);
	n += fwrite(series->coeff, sizeof(*series->coeff), sh_series_length(series->l_max, series->polar), file);

	return -(n != sh_series_length(series->l_max, series->polar) + 2);
}


/*
 * Read an sh_series object from a file.  Returns 0 on success, -1 on
 * failure (use ferror() or errno for details).
 *
 * This is a hack, don't use in production code.
 */


int sh_series_read(struct sh_series *series, FILE *file)
{
	unsigned int l_max;
	int polar;
	size_t n;

	n = fread(&l_max, sizeof(l_max), 1, file);
	n += fread(&polar, sizeof(polar), 1, file);
	sh_series_resize(series, l_max);
	sh_series_set_polar(series, polar);
	n += fread(series->coeff, sizeof(*series->coeff), sh_series_length(series->l_max, series->polar), file);

	return -(n != sh_series_length(series->l_max, series->polar) + 2);
}
