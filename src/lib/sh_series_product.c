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
 * the product there, and transform back to the frequency domain.
 *
 * The product plans can require a large quantity of memory.  They require
 * 20 bytes for every multiplication that must be performed to compute the
 * product, which scales as O(l^5) for non-azimuthally symmetric
 * multiplicands.  For l = 50, the plan for non-azimuthally symmetric
 * multiplicands requires over 6 GB of RAM.
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


double sh_series_wigner_3j(int ja, int jb, int jc, int ma, int mb, int mc)
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
 * Construct a product evaluation plan.  For internal use only.
 */


static int _product_plan(struct _sh_series_product_plan_op *microcode, int D_l_max, int D_polar, int A_l_max, int A_polar, int B_l_max, int B_polar)
{
	struct _sh_series_product_plan_op *op = microcode;
	int d_l, d_m, a_l, a_m, b_l;

	for(d_l = 0; d_l <= D_l_max; d_l++) {
		int d_m_max = D_polar || (A_polar && B_polar) ? 0 : d_l;
		for(d_m = -d_m_max; d_m <= d_m_max; d_m++) {
			int a_l_min = d_l > B_l_max ? d_l - B_l_max : 0;
			int a_l_max = A_l_max < d_l + B_l_max ? A_l_max : d_l + B_l_max;
			for(a_l = a_l_min; a_l <= a_l_max; a_l++) {
				int b_l_max = B_l_max < d_l + a_l ? B_l_max : d_l + a_l;
				int b_l_min = abs(d_l - a_l);
				b_l_min += (a_l + b_l_min + d_l) & 1;
				for(b_l = b_l_min; b_l <= b_l_max; b_l += 2) {
					int a_m_min = -(a_l < b_l - d_m ? a_l : b_l - d_m);
					int a_m_max = a_m_max = a_l < b_l + d_m ? a_l : b_l + d_m;
					if(A_polar) {
						if(a_m_min > 0 || a_m_max < 0)
							continue;
						a_m_min = a_m_max = 0;
					}
					for(a_m = a_m_min; a_m <= a_m_max; a_m++) {
						int b_m = d_m - a_m;
						if(B_polar && b_m != 0)
							continue;
						*op++ = (struct _sh_series_product_plan_op) {
							.dest_offset = sh_series_params_lmoffset(D_l_max, d_l, d_m),
							.a_offset = sh_series_params_lmoffset(A_l_max, a_l, a_m),
							.b_offset = sh_series_params_lmoffset(B_l_max, b_l, b_m),
							.factor = (d_m & 1 ? -1.0 : +1.0) * sqrt((2 * d_l + 1) * (2 * a_l + 1) * (2 * b_l + 1) / (4 * M_PI)) * sh_series_wigner_3j(a_l, b_l, d_l, 0, 0, 0) * sh_series_wigner_3j(a_l, b_l, d_l, a_m, b_m, -d_m)
						};
					}
				}
			}
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

	if(!new) {
		free(new);
		return NULL;
	}

	new->a_l_max = a->l_max;
	new->a_polar = a->polar;
	new->b_l_max = b->l_max;
	new->b_polar = b->polar;
	new->dest_l_max = dest_l_max;
	new->dest_polar = dest->polar;
	new->plan_length = 0;
	new->microcode = NULL;

	/* for small series, use frequency-domain algorithm, for large
	 * series use pixel-domain algorithm.  unfortunately "small" and
	 * "large", here, are defined by the stability of GSL's Wigner 3-j
	 * implementation, not by which approach is actually fastest.  we
	 * need to switch to the pixel domain implementation at smaller l's
	 * than those for which it is actually faster to do so */
	if(dest_l_max < 60) {
		struct _sh_series_product_plan_op *microcode;

		/* allocate worst-case size for microcode */
		new->microcode = malloc(sh_series_length(dest_l_max, dest->polar) * sh_series_length(a->l_max, a->polar) * sh_series_length(b->l_max, b->polar) * sizeof(*new->microcode));
		if(!new->microcode) {
			free(new);
			return NULL;
		}

		new->plan_length = _product_plan(new->microcode, dest_l_max, dest->polar, a->l_max, a->polar, b->l_max, b->polar);

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
	complex double *dst = dest->coeff;
	const complex double *alm = a->coeff;
	const complex double *blm = b->coeff;

	/* zero the destination */
	sh_series_zero(dest);

	/* execute the microcode */
	for(; op < last_op; op++)
		dst[op->dest_offset] += alm[op->a_offset] * blm[op->b_offset] * op->factor;

	return dest;
}


static struct sh_series *pixel_domain_product(struct sh_series *dest, const struct sh_series *a, const struct sh_series *b)
{
	unsigned l_max = a->l_max + b->l_max;	/* l_max of product */
	struct sh_series *wrkspc;
	complex double *a_mesh;
	complex double *b_mesh;
	int i;

	/* up-sample a and b inputs to l_max of product, and convert to
	 * pixel domain */

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

	/* compute the product in the pixel domain */

	/* FIXME:  don't assume we know the mesh size */
	for(i = 0; i < (int) (2 * (l_max + 1) * 2 * (l_max + 1)); i++)
		a_mesh[i] *= b_mesh[i];
	free(b_mesh);

	/* convert result to frequency domain */

	wrkspc = sh_series_new(l_max, a->polar && b->polar);
	if(!wrkspc || !sh_series_from_mesh(wrkspc, a_mesh)) {
		sh_series_free(wrkspc);
		free(a_mesh);
		return NULL;
	}
	free(a_mesh);

	/* copy into destination, possibly resizing to new l_max */

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
	if((plan->a_l_max != a->l_max) || (plan->b_l_max != b->l_max) || (plan->dest_l_max > dest->l_max) || (plan->a_polar != a->polar) || (plan->b_polar != b->polar) || (plan->dest_polar != dest->polar))
		return NULL;

	if(plan->microcode)
		return frequency_domain_product(dest, a, b, plan);
	return pixel_domain_product(dest, a, b);
}
