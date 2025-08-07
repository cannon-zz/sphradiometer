/*
 * Copyright (C) 2006--2009,2012,2019  Kipp C. Cannon
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
 * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 *
 * Do not, under any circumstances, edit this code.  Yes, this means you!
 * No, I'm really not kidding.  You cannot begin to imagine how much effort
 * I have gone to to check and cross-check the code in this file.  I have
 * never, ever, seen algorithms such as this that can be wrong in so many
 * different ways and still produce the correct answer.  The only way to
 * ensure this code works is to painstakingly construct test cases that
 * excercise every factor of I, and -1, and every conditional and ensure
 * that each is correct individually.  Please, I beg of you...
 *
 * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 */

/*
 * FIXME:  with the conventions adopted here, the special case of rotations
 * around the y axis produce real-valued D matrixes.  rotations around the
 * z axis, meanwhile, are trivially implemented as multiplications by
 * complex phase factors.  consider rewriting the general rotation case as
 * a three step algorithm consisting of a rotation about the z axis
 * (implemented as complex phase factors), followed by rotation about the y
 * axis (implemented with real-valued D matrixes), followed by a final
 * rotation about the z axis (again implemented as complex phase factors).
 * the D matrixes would require 1/2 the memory to store and altogether an
 * arbitrary rotation would require approximately 1/3 as many arithmetic
 * operations
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
 *                              Rotation Matrix
 *
 * ============================================================================
 */


/*
 * The elements of the 3x3 rotation matrix R are stored in the order
 *
 *	[ x' ]   [ R[0], R[1], R[2] ] [ x ]
 *	[ y' ] = [ R[3], R[4], R[5] ] [ y ]
 *	[ z' ]   [ R[6], R[7], R[8] ] [ z ]
 *
 * The element order is standard C order, and defining the rotation matrix
 * to be applied as indicated above allows the transformation to be
 * computed in a cache-friendly manner (all arrays are accessed
 * sequentially).
 */


/*
 * The matrix returned by euler_rotation_matrix() rotates by gamma about
 * the positive z axis, then by beta about the positive y axis, then by
 * alpha about the positive z axis.  euler_inv_rotation_matrix() returns
 * the inverse of the same matrix.
 *
 * Compared to wikipedia's article on Euler angles, this is the
 * "y-convention Geometrical definition", with our gamma = wikipedia's
 * alpha, our beta = wikipedia's beta, our alpha = wikipedia's gamma.
 *
 * https://en.wikipedia.org/wiki/Euler_angles
 */


double *euler_rotation_matrix(double gamma, double beta, double alpha)
{
	enum {
		x = 0,
		y = 1,
		z = 2
	};
	double *R = malloc(9 * sizeof(*R));
	double cos_alpha = cos(alpha);
	double sin_alpha = sin(alpha);
	double cos_beta = cos(beta);
	double sin_beta = sin(beta);
	double cos_gamma = cos(gamma);
	double sin_gamma = sin(gamma);

	if(!R)
		return NULL;

	/*
	 * using rotation matrixes from
	 *
	 * https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
	 *
	 * wxMaxima code to generate these expressions below:
	 *
	 * Rz(w) := matrix([cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]);
	 * Ry(w) := matrix([cos(w), 0, sin(w)], [0, 1, 0], [-sin(w), 0, cos(w)]);
	 *
	 * Rz(alpha).Ry(beta).Rz(gamma);
	 *
	 * The result is identical to wikipedia's right-handed, active
	 * rotation Z_alpha Y_beta Z_gamma result for proper Euler angles,
	 * and that is the convention adoped on wikipedia's page on the
	 * Wigner D-matrix, described there as the "z-y-z convention,
	 * right-handed frame, right-hand screw rule, active
	 * interpretation)".
	 */

	R[3 * x + x] =  cos_alpha * cos_beta * cos_gamma - sin_alpha * sin_gamma;
	R[3 * x + y] = -cos_alpha * cos_beta * sin_gamma - sin_alpha * cos_gamma;
	R[3 * x + z] =  cos_alpha * sin_beta;
	R[3 * y + x] =  sin_alpha * cos_beta * cos_gamma + cos_alpha * sin_gamma;
	R[3 * y + y] = -sin_alpha * cos_beta * sin_gamma + cos_alpha * cos_gamma;
	R[3 * y + z] =  sin_alpha * sin_beta;
	R[3 * z + x] = -sin_beta * cos_gamma;
	R[3 * z + y] =  sin_beta * sin_gamma;
	R[3 * z + z] =  cos_beta;

	return R;
}


double *euler_inv_rotation_matrix(double gamma, double beta, double alpha)
{
	return euler_rotation_matrix(-alpha, -beta, -gamma);
}


/*
 * The matrix returned by this function rotates the \hat{z} vector to point
 * in the direction given by the spherical polar co-ordinates theta and
 * phi.  (the rotation is by theta about the y axis, then by phi about the
 * z axis).
 */


double *sh_series_rot_matrix(double theta, double phi)
{
	return euler_rotation_matrix(0., theta, phi);
}


/*
 * The matrix returned by this function rotates the \hat{s} vector pointing
 * in the direction given by the spherical polar co-ordinates theta and phi
 * to point in the \hat{z} direction.  (the rotation is by
 * -phi about the z axis, then by -theta about the y axis).
 */


double *sh_series_invrot_matrix(double theta, double phi)
{
	return euler_inv_rotation_matrix(0., theta, phi);
}


/*
 * ============================================================================
 *
 *                Spherical Harmonic Series Rotation Matrices
 *
 * ============================================================================
 */


/*
 * The code in this section implements the algorithm of Choi, Ivanic,
 * Gordon, and Ruedenberg (Journal of Chemical Physics, 111 (19), Nov.
 * 1999, pg. 8825) for determining the rotation matrices for spherical
 * harmonics.  This is a recursive algorithm.  The matrix for l = 1 is
 * hard-coded, and the matrix for l > 1 is determined from the matrices for
 * l - 1 and 1 (l = 0 is a no-op).
 *
 * Note the following conventions.  Firstly, that the rotation matrix used
 * as input to the code here (see above) is the transpose of the definition
 * of the rotation matrix in Choi et al.  Secondly, in using a D matrix, a
 * spherical harmonic is transformed by summing over the m index,
 *
 * 	Y'_{l m'} = \sum_{m = -l}^{l} D^{l}_{m m'} Y_{l m},
 *
 * while a coefficient in an expansion is transformed by summing over the
 * m' index,
 *
 *	F'_{l m} = \sum_{m' = -l}^{l} D_{l}_{m m'} F_{l m'},
 *
 * where
 *
 * 	f(\theta, \phi) = \sum_{l, m} F_{l m} Y_{l m}(\theta, \phi).
 *
 * Finally, the DElement() macro (see below) defines the order in which the
 * D matrix components are stored in memory (the order is optimized for the
 * case of transforming coefficients in harmonic expansions of functions).
 *
 * This means that coefficients are transformed by summing across a row in
 * Choi et al.'s equations (5.4) and (5.5).
 *
 * Note also the error of an over-all sign in equation (5.5) of Choi et
 * al.
 */


/*
 * D matrix create, destroy and copy.  The pointer is shifted so that negative m
 * indeces can be used.
 */


static complex double *D_matrix_new(unsigned int l)
{
	complex double *D = malloc((2 * l + 1) * (2 * l + 1) * sizeof(*D));

	if(D)
		D += (2 * l + 1) * l + l;

	return D;
}


static void D_matrix_free(complex double *D, unsigned int l)
{
	if(D)
		free(D - (2 * l + 1) * l - l);
}


static complex double *D_matrix_copy(complex double *D, unsigned int l)
{
	complex double *new = malloc((2 * l + 1) * (2 * l + 1) * sizeof(*D));

	if(new) {
		memcpy(new, D - (2 * l + 1) * l - l, (2 * l + 1) * (2 * l + 1) * sizeof(*D));
		new += (2 * l + 1) * l + l;
	}

	return new;
}


/*
 * How to access an element in a D matrix.  The type cast to int is because
 * it is sometimes convenient to store l values as unsigned ints, but m's
 * are signed and can be negative, so we need to force the correct cast.
 */


#define DElement(matrix, l, m, m_prime) (matrix)[((int) (2 * (l) + 1)) * (m_prime) + (m)]


/*
 * Equations (5.3)--(5.5), the rotation matrix for l = 1.
 */


static complex double *D1(const double *R)
{
	enum {
		x = 0,
		y = 1,
		z = 2
	};
	complex double *D = D_matrix_new(1);

	if(!D)
		return NULL;

	DElement(D, 1, -1, -1) = (R[3 * y + y] + R[3 * x + x] - I * (R[3 * x + y] - R[3 * y + x])) / 2;
	DElement(D, 1, -1,  0) = (R[3 * z + x] - I * R[3 * z + y]) / sqrt(2);
	DElement(D, 1, -1, +1) = (R[3 * y + y] - R[3 * x + x] + I * (R[3 * x + y] + R[3 * y + x])) / 2;
	DElement(D, 1,  0, -1) = (R[3 * x + z] + I * R[3 * y + z]) / sqrt(2);
	DElement(D, 1,  0,  0) = R[3 * z + z];
	DElement(D, 1,  0, +1) = (-R[3 * x + z] + I * R[3 * y + z]) / sqrt(2);
	DElement(D, 1, +1, -1) = (R[3 * y + y] - R[3 * x + x] - I * (R[3 * x + y] + R[3 * y + x])) / 2;
	DElement(D, 1, +1,  0) = (-R[3 * z + x] - I * R[3 * z + y]) / sqrt(2);
	DElement(D, 1, +1, +1) = (R[3 * y + y] + R[3 * x + x] - I * (R[3 * y + x] - R[3 * x + y])) / 2;

	return D;
}


/*
 * The recursion given in (6.1), (6.9), and (6.14).
 */


static double a(int l, int m, int m_prime)
{
	/* (6.2) */
	return sqrt((l + m) * (l - m) / (double) ((l + m_prime) * (l - m_prime)));
}


static double b(int l, int m, int m_prime)
{
	/* (6.3) */
	return sqrt((l + m) * (l + m - 1) / (double) (2 * (l + m_prime) * (l - m_prime)));
}


static double c(int l, int m, int m_prime)
{
	/* (6.10) */
	return sqrt(2 * (l + m) * (l - m) / (double) ((l + m_prime) * (l + m_prime - 1)));
}


static double d(int l, int m, int m_prime)
{
	/* (6.11) */
	return sqrt((l + m) * (l + m - 1) / (double) ((l + m_prime) * (l + m_prime - 1)));
}


static complex double *DL(int l, const complex double *DL_minus_1, const complex double *D1)
{
	complex double *D = D_matrix_new(l);
	int m, m_prime;

	if(!D)
		return NULL;

	for(m = -l; m <= l; m++) {
		/* (6.1), (6.4), and (6.5) */
		for(m_prime = -l + 1; m_prime <= l - 1; m_prime++)
			DElement(D, l, m, m_prime) = (abs(m) == l ? 0 : a(l, m, m_prime) * DElement(D1, 1, 0, 0) * DElement(DL_minus_1, l - 1, m, m_prime)) + (m <= -l + 1 ? 0 : b(l, m, m_prime) * DElement(D1, 1, 1, 0) * DElement(DL_minus_1, l - 1, m - 1, m_prime)) + (-m <= -l + 1 ? 0 : b(l, -m, m_prime) * DElement(D1, 1, -1, 0) * DElement(DL_minus_1, l - 1, m + 1, m_prime));

		/* (6.9), (6.12), and (6.13) */
		m_prime = -l;
		DElement(D, l, m, m_prime) = (abs(m) == l ? 0 : c(l, m, -m_prime) * DElement(D1, 1, 0, -1) * DElement(DL_minus_1, l - 1, m, m_prime + 1)) + (m <= -l + 1 ? 0 : d(l, m, -m_prime) * DElement(D1, 1, 1, -1) * DElement(DL_minus_1, l - 1, m - 1, m_prime + 1)) + (-m <= -l + 1 ? 0 : d(l, -m, -m_prime) * DElement(D1, 1, -1, -1) * DElement(DL_minus_1, l - 1, m + 1, m_prime + 1));

		/* (6.14), (6.12), and (6.13) */
		m_prime = l;
		DElement(D, l, m, m_prime) = (abs(m) == l ? 0 : c(l, m, m_prime) * DElement(D1, 1, 0, 1) * DElement(DL_minus_1, l - 1, m, m_prime - 1)) + (m <= -l + 1 ? 0 : d(l, m, m_prime) * DElement(D1, 1, 1, 1) * DElement(DL_minus_1, l - 1, m - 1, m_prime - 1)) + (-m <= -l + 1 ? 0 : d(l, -m, m_prime) * DElement(D1, 1, -1, 1) * DElement(DL_minus_1, l - 1, m + 1, m_prime - 1));
	}

	return D;
}


/*
 * ============================================================================
 *
 *                  Spherical Harmonic Series Rotation Plan
 *
 * ============================================================================
 */


/*
 * Construct a rotation plan.
 */


struct sh_series_rotation_plan *sh_series_rotation_plan_new(const struct sh_series *series, const double *R)
{
	struct sh_series_rotation_plan *plan = malloc(sizeof(*plan));
	complex double **D = malloc(series->l_max * sizeof(*D));
	unsigned int l = 1;

	/* if l_max = 0 then allow D to be NULL '*/
	if(!plan || (series->l_max > 0 && !D))
		goto error;

	/* adjust so D[1] is the first matrix */
	D--;

	/* apply recursive procedure to compute D matrices */
	if(series->l_max >= 1) {
		D[l] = D1(R);
		if(!D[l])
			goto error;
	}
	for(l++; l <= series->l_max; l++) {
		D[l] = DL(l, D[l - 1], D[1]);
		if(!D[l])
			goto error;
	}

	plan->l_max = series->l_max;
	plan->D = D;

	return plan;

error:
	while(--l >= 1)
		D_matrix_free(D[l], l);
	free(D + 1);
	free(plan);
	return NULL;
}


/*
 * Free a rotation plan.
 */


void sh_series_rotation_plan_free(struct sh_series_rotation_plan *plan)
{
	unsigned int l;

	if(plan) {
		for(l = 1; l <= plan->l_max; l++)
			D_matrix_free(plan->D[l], l);
		free(plan->D + 1);
	}
	free(plan);
}


/*
 * Copy a rotation plan.
 */


struct sh_series_rotation_plan *sh_series_rotation_plan_copy(const struct sh_series_rotation_plan *plan)
{
	unsigned l;
	struct sh_series_rotation_plan *new = malloc(sizeof(*new));

	*new = *plan;
	new->D = malloc(plan->l_max * sizeof(*new->D));
	/* l starts from 1. fit a start position of D with the index l. */
	new->D--;
	for(l = 1; l <= plan->l_max; l++)
		new->D[l] = D_matrix_copy(plan->D[l], l);

	return new;
}


/*
 * Change the l_max of a rotation plan.  The conversion is done in-place.
 * The return value is NULL on failure, or the input plan on success.  On
 * failure, the input plan is left unchanged.
 *
 * NOTE:  for now only reduction of l_max is supported.
 */


struct sh_series_rotation_plan *sh_series_rotation_plan_set_l(struct sh_series_rotation_plan *plan, unsigned l_max)
{
	unsigned l;
	complex double **D, **Dorig;

	/* error? */
	/* FIXME:  generalize to allow l_max to increase */
	if(l_max > plan->l_max)
		return NULL;
	/* no-op? */
	if(l_max == plan->l_max)
		return plan;

	/* save the contents of the array of D matrix addresses so we can
	 * free the matrices if the realloc() succeeds.  we've promised to
	 * leave the plan unchanged if the realloc() fails so we can't free
	 * the matrices until we see that it succeeds */
	Dorig = malloc((plan->l_max - l_max) * sizeof(*Dorig));
	if(!Dorig)
		return NULL;
	memcpy(Dorig, &plan->D[l_max + 1], (plan->l_max - l_max) * sizeof(*Dorig));

	/* shorten the array of D matrix adddresses.  note that plan->D
	 * points to the memory address immediate preceding the actual
	 * allocated address because we don't need D[0] and don't want to
	 * have to screw with adding an offset all the time */
	D = realloc(plan->D + 1, l_max * sizeof(*D));
	if(!D) {
		free(Dorig);
		return NULL;
	}
	plan->D = D - 1;

	/* free the unused matrices */
	for(l = l_max + 1; l <= plan->l_max; l++)
		D_matrix_free(Dorig[l - l_max - 1], l);
	free(Dorig);

	plan->l_max = l_max;
	return plan;
}


/*
 * Access a Wigner D matrix element in an sh_series_rotation_plan.  Returns
 * the matrix element or complex NaN if the indexes do not correspond to a
 * valid element.
 */


complex double sh_series_rotation_plan_wigner_D(const struct sh_series_rotation_plan *plan, unsigned l, int m, int m_prime)
{
	if(l > plan->l_max || (unsigned) abs(m) > l || (unsigned) abs(m_prime) > l)
		return nan("") + I * nan("");

	return l == 0 ? 1.0 : DElement(plan->D[l], l, m, m_prime);
}


/*
 * ============================================================================
 *
 *                    Spherical Harmonic Series Rotations
 *
 * ============================================================================
 */


/*
 * Apply a rotation plan.
 */


static struct sh_series *_sh_series_rotate(struct sh_series *result, const struct sh_series *series, const struct sh_series_rotation_plan *plan)
{
	unsigned int l;

	/* parallelizing the outer loop via omp results in a 50% slow down
	 * of the test suite, with 4x CPU use. */
	for(l = 1; l <= series->l_max; l++) {
		complex double *D = plan->D[l];
		for(int m_prime = -(int) l; m_prime <= (int) l; m_prime++) {
			complex double c = 0.0;
			for(int m = -(int) l; m <= (int) l; m++)
				c += DElement(D, l, m, m_prime) * sh_series_get(series, l, m);
			sh_series_set(result, l, m_prime, c);
		}
	}

	return result;
}


static struct sh_series *_sh_series_rotate_polar(struct sh_series *result, const struct sh_series *series, const struct sh_series_rotation_plan *plan)
{
	unsigned int l;

	for(l = 1; l <= series->l_max; l++) {
		complex double *D = plan->D[l];
		for(int m_prime = -(int) l; m_prime <= (int) l; m_prime++)
			sh_series_set(result, l, m_prime, DElement(D, l, 0, m_prime) * sh_series_get(series, l, 0));
	}

	return result;
}


struct sh_series *sh_series_rotate(struct sh_series *result, const struct sh_series *series, const struct sh_series_rotation_plan *plan)
{
	if((result->l_max != series->l_max) || result->polar || (result == series))
		return NULL;

	sh_series_set(result, 0, 0, sh_series_get(series, 0, 0));
	if(series->polar)
		return _sh_series_rotate_polar(result, series, plan);
	return _sh_series_rotate(result, series, plan);
}


/*
 * Rotate about positive z-axis by angle phi.  This is a special case that
 * can be performed quickly, and without first computing any D matrices.
 * This transformation can also be done in place (the input and output
 * pointers can point to the same sh_series object).
 */


struct sh_series *sh_series_rotate_z(struct sh_series *result, const struct sh_series *series, double phi)
{
	if(result->l_max != series->l_max)
		return NULL;

	if(series->polar) {
		/* no-op */
		if(result != series)
			/* copy coefficients and zero any unused */
			sh_series_assign(result, series);
	} else {
		int m;

		if(result->polar)
			return NULL;

		for(m = -(int) series->l_max; m <= (int) series->l_max; m++) {
			const complex double factor = cexp(-I * m * phi);
			unsigned int l;
			for(l = abs(m); l <= series->l_max; l++)
				sh_series_set(result, l, m, sh_series_get(series, l, m) * factor);
		}
	}

	return result;
}
