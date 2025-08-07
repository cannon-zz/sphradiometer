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


#ifndef __RADIOMETER_SH_SERIES_H__
#define __RADIOMETER_SH_SERIES_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <complex.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#define complex _Complex
#endif


/*
 * ============================================================================
 *
 *                                 Data Types
 *
 * ============================================================================
 */


/*
 * An expansion of a function on the sphere, expressed as the coefficients
 * in a series of spherical harmonics upto and including order l = l_max.
 * The ordering of the coefficients in the array is as follows (expressed
 * as l,m pairs, read left-to-right, top-to-bottom):
 *
 * (0, 0), (1, 0), (2, 0), (3, 0), ... , (l_max, 0),
 *         (1,+1), (2,+1), (3,+1), ... , (l_max,+1),
 *                 (2,+2), (3,+2), ... , (l_max,+2),
 *                         (3,+3), ... , (l_max,+3),
 *                                        ...
 *                                       (l_max,+l_max),
 *                                       (l_max,-l_max),
 *                                        ...
 *                         (3,-3), ... , (l_max,-3),
 *                 (2,-2), (3,-2), ... , (l_max,-2),
 *         (1,-1), (2,-1), (3,-1), ... , (l_max,-1)
 *
 * In other words coefficients are arranged in order of increasing l for
 * fixed m.  This ordering allows for efficient summing of the series,
 * since there exist efficient recursive algorithms for evaluating
 * spherical harmonics over a range of m at fixed l.  Also, since
 * real-valued functions can be described by specifying only the m >= 0
 * coefficients, they can be represented using only the first half of the
 * coefficient array.  As well, if the "polar" flag is true --- the
 * function has azimuthal symmetry --- then only the m = 0 coefficients are
 * present.  But in all cases the coefficients that are present are located
 * at the same positions in memory allowing efficient copying and looping.
 *
 * This coefficient ordering is the same as used in S2Kit.
 */


struct sh_series {
	unsigned int l_max;
	int polar;
	double complex *coeff;
};


/*
 * Series product plan.  An opaque data type used to store pre-computed
 * coefficients for use in computing the product of two spherical harmonic
 * series as a spherical harmonic series.
 */


struct sh_series_product_plan {
	unsigned int a_l_max, b_l_max, dest_l_max;
	int a_polar, b_polar, dest_polar;
	int plan_length;
	struct _sh_series_product_plan_op {
		double factor;
		int dest_offset;
		int a_offset;
		int b_offset;
	} *microcode;
};


/*
 * Series rotation plan.  An opaque data type used to store the D matrix
 * elements for rotating a function expanded in a spherical harmonic
 * series.
 */


struct sh_series_rotation_plan {
	unsigned int l_max;
	double complex **D;
};


/*
 * sh_series_array object --- an array of sh_series objects sharing a
 * coefficient buffer.  The coefficients for each sh_series object are
 * stored contiguously in the buffer, with each sh_series object's
 * coefficients coming immediately after the coefficients for the object
 * preceding it in the array.
 */


struct sh_series_array {
	unsigned int l_max;
	int polar;
	int n;
	struct sh_series *series;
	double complex *coeff;
	int stride;
};


/*
 * opaque interpolator object
 */


struct sh_series_eval_interp;


/*
 * ============================================================================
 *
 *                                   Macros
 *
 * ============================================================================
 */


/*
 * Return the offset for the coefficient (l,m) = (abs(m), m).  Recalling
 * that the coefficients are grouped by m, the offset returned by this
 * function is the start of the coefficients with the given m (see above).
 */


#ifndef SWIG
static size_t sh_series_moffset(unsigned int, int) __attribute__ ((unused));
#endif
static size_t sh_series_moffset(unsigned int l_max, int m)
{
	/* C allows compilers to perform any re-ordering of algebraic
	 * operations permited by the normal rules of algebra.  here we
	 * rely on a specific evaluation order to get the integer
	 * arithmetic to yield the correct result.  to force the
	 * multiplication to be performed before the divide-by-2 that
	 * follows it, it has to be done this way, storing the result of
	 * the multiplication in an intermediate variable */
	const int x = m * (m - 1);
	l_max += 1;
	if(m >= 0)
		return l_max * m - x / 2;
	return l_max * (l_max + m) + x / 2;

}


/*
 * Return the offset of the (l,m)-th coefficient.
 */


#ifndef SWIG
static size_t sh_series_params_lmoffset(unsigned int, unsigned int, int) __attribute__ ((unused));
#endif
static size_t sh_series_params_lmoffset(unsigned int l_max, unsigned int l, int m)
{
	return sh_series_moffset(l_max, m) + l - abs(m);
}


#ifndef SWIG
static size_t sh_series_lmoffset(const struct sh_series *, unsigned int, int) __attribute__ ((unused));
#endif
static size_t sh_series_lmoffset(const struct sh_series *series, unsigned int l, int m)
{
	return sh_series_params_lmoffset(series->l_max, l, m);
}


/*
 * Return the total number of coefficients in the expansion of order l_max.
 */


#ifndef SWIG
static size_t sh_series_length(unsigned int, int) __attribute__ ((unused));
#endif
static size_t sh_series_length(unsigned int l_max, int polar)
{
	l_max += 1;
	if(polar)
		return l_max;
	return l_max * l_max;
}


/*
 * ============================================================================
 *
 *                               Misc Utilities
 *
 * ============================================================================
 */


/* wrapper of GSL's Wigner 3-j function */
double sh_series_wigner_3j(int, int, int, int, int, int);


/*
 * ============================================================================
 *
 *                   Spherical Harmonic Function Prototypes
 *
 * ============================================================================
 */


double complex sh_series_Y(unsigned int, int, double, double);
double complex sh_series_Yconj(unsigned int, int, double, double);
double complex *sh_series_Y_array(double complex *, unsigned int, int, double, double);
double complex *sh_series_Yconj_array(double complex *, unsigned int, int, double, double);


/*
 * ============================================================================
 *
 *                            sh_series Prototypes
 *
 * ============================================================================
 */


/*
 * Basics
 */


struct sh_series *sh_series_new(unsigned int, int);
struct sh_series *sh_series_new_zero(unsigned int, int);
struct sh_series *sh_series_copy(const struct sh_series *);
struct sh_series *sh_series_assign(struct sh_series *, const struct sh_series *);
struct sh_series *sh_series_resize(struct sh_series *, unsigned int);
struct sh_series *sh_series_zero(struct sh_series *);
struct sh_series *sh_series_set_polar(struct sh_series *, int);
void sh_series_free(struct sh_series *);


double complex sh_series_get(const struct sh_series *, unsigned int, int);
double complex sh_series_set(struct sh_series *, unsigned int, int, double complex);


/*
 * Harmonic domain <--> spatial domain transformations
 */


double complex sh_series_eval(const struct sh_series *, double, double);
struct sh_series_eval_interp *sh_series_eval_interp_new(const struct sh_series *);
void sh_series_eval_interp_free(struct sh_series_eval_interp *);
double complex sh_series_eval_interp(const struct sh_series_eval_interp *, double, double);


complex double *sh_series_mesh_new(unsigned int, int *, int *, double **, double **);
double *sh_series_real_mesh_new(unsigned int, int *, int *, double **, double **);
struct sh_series *sh_series_from_mesh(struct sh_series *, double complex *);
struct sh_series *sh_series_from_realmesh(struct sh_series *, double *);
double complex *sh_series_to_mesh(const struct sh_series *);


struct sh_series *sh_series_from_func(struct sh_series *, double complex (*)(double, double, void *), void *);
struct sh_series *sh_series_from_realfunc(struct sh_series *, double (*)(double, double, void *), void *);

struct sh_series *sh_series_impulse(unsigned int, double, double);


/*
 * Arithmetic
 */


struct sh_series *sh_series_add(struct sh_series *, double complex, const struct sh_series *);
struct sh_series *sh_series_scale(struct sh_series *, double complex);
struct sh_series *sh_series_conj(struct sh_series *);
struct sh_series *sh_series_real(struct sh_series *);
struct sh_series *sh_series_imag(struct sh_series *);
struct sh_series *sh_series_clip(struct sh_series *, double);
struct sh_series *sh_series_scale_power_l(struct sh_series *, double, double);
double complex sh_series_dot(const struct sh_series *, const struct sh_series *);

struct sh_series_product_plan *sh_series_product_plan_new(const struct sh_series *, const struct sh_series *, const struct sh_series *);
void sh_series_product_plan_free(struct sh_series_product_plan *);
struct sh_series *sh_series_product(struct sh_series *, const struct sh_series *, const struct sh_series *, const struct sh_series_product_plan *);
int sh_series_is_nan(const struct sh_series *);
int sh_series_is_inf(const struct sh_series *);


/*
 * Rotation
 */


double *euler_rotation_matrix(double, double, double);
double *euler_inv_rotation_matrix(double, double, double);
double *sh_series_rot_matrix(double, double);
double *sh_series_invrot_matrix(double, double);
struct sh_series_rotation_plan *sh_series_rotation_plan_new(const struct sh_series *, const double *);
void sh_series_rotation_plan_free(struct sh_series_rotation_plan *);
struct sh_series_rotation_plan *sh_series_rotation_plan_copy(const struct sh_series_rotation_plan *);
struct sh_series_rotation_plan *sh_series_rotation_plan_set_l(struct sh_series_rotation_plan *, unsigned);
double complex sh_series_rotation_plan_wigner_D(const struct sh_series_rotation_plan *, unsigned, int, int);
struct sh_series *sh_series_rotate(struct sh_series *, const struct sh_series *, const struct sh_series_rotation_plan *);
struct sh_series *sh_series_rotate_z(struct sh_series *, const struct sh_series *, double);


/*
 * Differentiation and integration
 */


struct sh_series *sh_series_times_costheta(const struct sh_series *);
struct sh_series *sh_series_d_by_dphi(struct sh_series *);
struct sh_series *sh_series_sintheta_d_by_dtheta(const struct sh_series *);
const struct sh_series *sh_series_sintheta_grad(const struct sh_series *, struct sh_series **, struct sh_series **);
struct sh_series *sh_series_sintheta_div(const struct sh_series *, const struct sh_series *);
struct sh_series *sh_series_sintheta_curl(const struct sh_series *, const struct sh_series *);
struct sh_series *sh_series_laplacian(struct sh_series *);
struct sh_series *sh_series_invlaplacian(struct sh_series *);


/*
 * I/O
 */


void sh_series_print(FILE *, const struct sh_series *);
int sh_series_write_healpix_map(const struct sh_series *, const char *);
int sh_series_write_healpix_alm(const struct sh_series *, const char *);
struct sh_series *sh_series_read_healpix_alm(const char *);


/*
 * ============================================================================
 *
 *                         sh_series_array Prototypes
 *
 * ============================================================================
 */


/*
 * Basics
 */


struct sh_series_array *sh_series_array_new(int, unsigned int, int);
void sh_series_array_free(struct sh_series_array *);
struct sh_series_array *sh_series_array_set_len(struct sh_series_array *, int);
struct sh_series_array *sh_series_array_set_polar(struct sh_series_array *, int);
struct sh_series_array *sh_series_array_set_l(struct sh_series_array *, unsigned int);
struct sh_series_array *sh_series_array_copy(const struct sh_series_array *);
struct sh_series_array *sh_series_array_assign(struct sh_series_array *, const struct sh_series_array *);
struct sh_series_array *sh_series_array_scale(struct sh_series_array *, double complex);


/*
 * Arithmetic
 */


struct sh_series *sh_series_array_dot(struct sh_series *, const struct sh_series_array *, const double *);
struct sh_series *sh_series_array_dotc(struct sh_series *, const struct sh_series_array *, const double complex *);


struct sh_series_array *sh_series_array_window(struct sh_series_array *, const double *);
struct sh_series_array *sh_series_array_windowc(struct sh_series_array *, const double complex *);


/*
 * Fourier transforms
 */


struct sh_series_array *sh_series_array_forward_fft(struct sh_series_array *);
struct sh_series_array *sh_series_array_reverse_fft(struct sh_series_array *);


#ifdef __cplusplus
#undef complex
}
#endif


#endif	/* __RADIOMETER_SH_SERIES_H__ */
