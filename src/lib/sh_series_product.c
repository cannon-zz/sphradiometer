/*
 * Copyright (C) 2006--2009,2019  Kipp C. Cannon
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


/*
 * Notes.  Computing the product of two functions known only by the
 * coefficients of their expansion in spherical harmonics is very costly.
 * The naive implementation is O(l_max^6).  By taking advantage of known
 * 0's of the Wigner 3-j symbols, it is possible to significantly reduce
 * the number of operations.  If the series are known to contain only m=0
 * terms, that is they are azimuthally-symmetric, then even greater
 * reduction of operation count is possible, but even in this extreme case
 * the algorithm is still O(l_max^3).
 *
 * There exist fast spherical harmonic transforms (FSHTs) in analogy to the
 * fast Fourier transform.  The leading co-efficients are very large.
 * Nevertheless, they exhibit asymptotic complexity (operation counts) of
 * O(l_max^2 log^2 l_max), which means that for large l it is
 * computationally advantageous to transform to the spatial domain, compute
 * the product there, and transform back to the frequency domain.  This is
 * not done here.
 *
 * The product plans can require a large quantity of memory.  They require
 * 20 bytes for every multiplication that must be performed to compute the
 * product, which scales as O(l^5) for non-azimuthally symmetric
 * multiplicands.  For l = 50, the plan for non-azimuthally symmetric
 * multiplicands requires over 6 GB of RAM.  But multiplications of that
 * size are absurd to perform in the harmonic domain;  the time should be
 * spent on incorporating an FSHT for those problems.
 */


#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                 Utilities
 *
 * ============================================================================
 */


/*
 * Wrapper for GSL's Wigner 3-j symbol function.  The 3-j symbols are
 * defined for integer and 1/2 integer arguments.  GSL uses integer
 * arguments, so to handle 1/2 integer values it uses the convention of
 * expecting its inputs to be doubled.
 */


static double wigner_3j(int ja, int jb, int jc, int ma, int mb, int mc)
{
	return gsl_sf_coupling_3j(2 * ja, 2 * jb, 2 * jc, 2 * ma, 2 * mb, 2 * mc);
}


/*
 * Compare two product plan microcode operations by the array offsets they
 * reference.  Used to sort a list of operations so as to access memory in
 * a cache-friendly manner.
 */


static int op_cmp(const void *arg1, const void *arg2)
{
	const struct _sh_series_product_plan_op *op1 = arg1;
	const struct _sh_series_product_plan_op *op2 = arg2;
	int result;

	result = op1->dest_offset - op2->dest_offset;
	if(!result) {
		result = op1->a_offset - op2->a_offset;
		if(!result)
			result = op1->b_offset - op2->b_offset;
	}
	return result;
}


/*
 * ============================================================================
 *
 *                 sh_series Object Product Evaluation Plans
 *
 * ============================================================================
 */


/*
 * Construct a product evaluation plan for non-azimuthally symmetric
 * multiplicands.  For internal use only.
 */


static int _product_plan(struct _sh_series_product_plan_op *microcode, int dest_l_max, int a_l_max, int b_l_max)
{
	struct _sh_series_product_plan_op *op = microcode;
	int p, q, k, m, l;

	for(p = 0; p <= dest_l_max; p++)
		for(q = -p; q <= p; q++) {
			const int kmax = a_l_max < p + b_l_max ? a_l_max : p + b_l_max;
			for(k = p > b_l_max ? p - b_l_max : 0; k <= kmax; k++) {
				const int mmax = b_l_max < p + k ? b_l_max : p + k;
				int mmin = abs(p - k);
				mmin += (k + mmin + p) & 1;
				for(m = mmin; m <= mmax; m += 2) {
					const int lmax = k < q + m ? k : q + m;
					for(l = -k > q - m ? -k : q - m; l <= lmax; l++)
						*op++ = (struct _sh_series_product_plan_op) {
							.dest_offset = sh_series_params_lmoffset(dest_l_max, p, q),
							.a_offset = sh_series_params_lmoffset(a_l_max, k, l),
							.b_offset = sh_series_params_lmoffset(b_l_max, m, q - l),
							.factor = (q & 1 ? -1.0 : +1.0) * sqrt((2 * p + 1) * (2 * k + 1) * (2 * m + 1) / (4 * M_PI)) * wigner_3j(k, m, p, 0, 0, 0) * wigner_3j(k, m, p, l, q - l, -q)
						};
				}
			}
		}

	return op - microcode;
}


/*
 * Construct a product evaluation plan for azimuthally symmetric
 * multiplicands.  For internal use only.
 */


static int _product_plan_polar(struct _sh_series_product_plan_op *microcode, int dest_l_max, int a_l_max, int b_l_max)
{
	struct _sh_series_product_plan_op *op = microcode;
	int p, k, m;

	for(p = 0; p <= dest_l_max; p++) {
		const int kmax = a_l_max < p + b_l_max ? a_l_max : p + b_l_max;
		for(k = p > b_l_max ? p - b_l_max : 0; k <= kmax; k++) {
			const int mmax = b_l_max < p + k ? b_l_max : p + k;
			int mmin = abs(p - k);
			mmin += (k + mmin + p) & 1;
			for(m = mmin; m <= mmax; m += 2)
				*op++ = (struct _sh_series_product_plan_op) {
					.dest_offset = sh_series_params_lmoffset(dest_l_max, p, 0),
					.a_offset = sh_series_params_lmoffset(a_l_max, k, 0),
					.b_offset = sh_series_params_lmoffset(b_l_max, m, 0),
					.factor = sqrt((2 * p + 1) * (2 * k + 1) * (2 * m + 1) / (4 * M_PI)) * wigner_3j(k, m, p, 0, 0, 0) * wigner_3j(k, m, p, 0, 0, 0)
				};
		}
	}

	return op - microcode;
}


/*
 * Construct a new product evaluation plan.  The "a" and "b" multiplicands
 * must have the same symmetry setting (azimuthally or non-azimuthally
 * symmetric), and if they are not azimuthally symmetric then neither must
 * be the product destination (although if the multiplicands are
 * azimuthally symmetric, the destination need not be).
 */


struct sh_series_product_plan *sh_series_product_plan_new(const struct sh_series *dest, const struct sh_series *a, const struct sh_series *b)
{
	struct sh_series_product_plan *new = malloc(sizeof(*new));
	const int dest_l_max = dest->l_max <= a->l_max + b->l_max ? dest->l_max : a->l_max + b->l_max;
	const int polar = a->polar && b->polar;
	struct _sh_series_product_plan_op *microcode;

	if(!new || (a->polar != b->polar) || (dest->polar && !polar)) {
		free(new);
		return NULL;
	}

	if((a->l_max > b->l_max ? a->l_max : b->l_max) < 40) {
		/* for small series, use frequency-domain algorithm */
		/* allocate worst-case size for microcode */
		microcode =  malloc(sh_series_length(dest_l_max, polar) * sh_series_length(a->l_max, polar) * sh_series_length(b->l_max, polar) * sizeof(*microcode));
		if(!microcode) {
			free(new);
			return NULL;
		}
	} else {
		/* for large series use pixel-domain algorithm */
		microcode = NULL;
	}

	new->a_l_max = a->l_max;
	new->b_l_max = b->l_max;
	new->dest_l_max = dest_l_max;
	new->polar = polar;
	new->plan_length = 0;
	new->microcode = microcode;

	if(microcode) {
		if(new->polar)
			new->plan_length = _product_plan_polar(microcode, dest_l_max, a->l_max, b->l_max);
		else
			new->plan_length = _product_plan(microcode, dest_l_max, a->l_max, b->l_max);

		/* shrink microcode to actual size */
		microcode = realloc(new->microcode, new->plan_length * sizeof(*new->microcode));
		if(!microcode) {
			free(new->microcode);
			free(new);
			return NULL;
		}
		new->microcode = microcode;

		/* sort the operations to optimize cache usage */
		qsort(new->microcode, new->plan_length, sizeof(*new->microcode), op_cmp);
	}

	return new;
}


/*
 * Destroy a product plan.
 */


void sh_series_product_plan_free(struct sh_series_product_plan *plan)
{
	if(plan)
		free(plan->microcode);
	free(plan);
}


/*
 * ============================================================================
 *
 *                             Product Evaluation
 *
 * ============================================================================
 */


static struct sh_series *frequency_domain_product(struct sh_series *dest, const struct sh_series *a, const struct sh_series *b, const struct sh_series_product_plan *plan)
{
	const struct _sh_series_product_plan_op *op = plan->microcode;
	const struct _sh_series_product_plan_op *last_op = plan->microcode + plan->plan_length;

	/* zero the destination */
	sh_series_zero(dest);

	/* execute the microcode */
	for(; op < last_op; op++)
		dest->coeff[op->dest_offset] += a->coeff[op->a_offset] * b->coeff[op->b_offset] * op->factor;

	return dest;
}


static struct sh_series *pixel_domain_product(struct sh_series *dest, const struct sh_series *a, const struct sh_series *b)
{
	unsigned l_max = a->l_max + b->l_max;
	struct sh_series *wrkspc;
	complex double *a_mesh;
	complex double *b_mesh;
	int i;

	wrkspc = sh_series_copy(a);
	if(!wrkspc || !sh_series_resize(wrkspc, l_max)) {
		sh_series_free(wrkspc);
		return NULL;
	}
	a_mesh = sh_series_to_mesh(wrkspc);
	sh_series_free(wrkspc);

	wrkspc = sh_series_copy(b);
	if(!wrkspc || !sh_series_resize(wrkspc, l_max)) {
		sh_series_free(wrkspc);
		free(a_mesh);
		return NULL;
	}
	b_mesh = sh_series_to_mesh(wrkspc);
	sh_series_free(wrkspc);

	if(!a_mesh || !b_mesh) {
		free(a_mesh);
		free(b_mesh);
		return NULL;
	}

	/* FIXME:  don't assume we know the mesh size */
	for(i = 0; i < (int) (2 * (l_max + 1) * 2 * (l_max + 1)); i++)
		a_mesh[i] *= b_mesh[i];
	free(b_mesh);

	wrkspc = sh_series_new(l_max, a->polar && b->polar);
	if(!wrkspc || !sh_series_from_mesh(wrkspc, a_mesh)) {
		sh_series_free(wrkspc);
		free(a_mesh);
		return NULL;
	}
	free(a_mesh);
	if(!sh_series_resize(wrkspc, dest->l_max) || !sh_series_assign(dest, wrkspc)) {
		sh_series_free(wrkspc);
		return NULL;
	}
	sh_series_free(wrkspc);

	return dest;
}


/*
 * Compute the product of two sh_series objects.  a and b can be the same
 * series, but dest must not point to either of them.
 */


struct sh_series *sh_series_product(struct sh_series *dest, const struct sh_series *a, const struct sh_series *b, const struct sh_series_product_plan *plan)
{
	/* check that the plan is appropriate */
	/* FIXME:  this isn't really needed if we're using the pixel-domain
	 * implementation */
	if((plan->a_l_max != a->l_max) || (plan->b_l_max != b->l_max) || (plan->dest_l_max > dest->l_max) || (plan->polar != a->polar) || (plan->polar != b->polar) || (!plan->polar && dest->polar))
		return NULL;

	if(plan->microcode)
		return frequency_domain_product(dest, a, b, plan);
	return pixel_domain_product(dest, a, b);
}

