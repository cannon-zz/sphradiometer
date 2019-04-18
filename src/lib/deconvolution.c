/*
 * Copyright (C) 2006--2009,2012,2019  Kipp C. Cannon
 * Copyright (C) 2010  Dan Stratman
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
#include <lapacke.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/deconvolution.h>


/*
 * ============================================================================
 *
 *                Local Wrapper for LAPACK's SVD Factorization
 *
 * ============================================================================
 */


static complex double *transpose(int rows, int cols, complex double **matrix)
{
	complex double *new = malloc(rows * cols * sizeof(*new));
	int i, j;

	if(!new)
		return NULL;

	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			new[j * rows + i] = (*matrix)[i * cols + j];

	free(*matrix);
	*matrix = new;

        return new;
}


struct SVD {
	int m;
	int n;
	complex double *A;
	complex double *A_new;
	complex double *U;
	complex double *VT;
	double *S;
};


static void SVD_free(struct SVD *singular_value_decomp)
{
	if(singular_value_decomp) {
		free(singular_value_decomp->A);
		free(singular_value_decomp->A_new);
		free(singular_value_decomp->U);
		free(singular_value_decomp->VT);
		free(singular_value_decomp->S);
	}
	free(singular_value_decomp);
}


static struct SVD *SVD_new(int m, int n, const complex double *A)
{
	struct SVD *new = malloc(sizeof(*new));

	int max = m > n ? m : n;
	int min = m < n ? m : n;

	char job = 'S';
	int lda = max;
	int ldu = m;
	int ldvt = min;
	int lwork = -1;
	int lrwork = 6 * min * min + 7 * min;
	int info = 0;

	if(!new)
		return NULL;

	new->A = malloc(m * n * sizeof(*new->A));
	new->A_new = malloc(m * n * sizeof(*new->A_new));
	new->S = calloc(min, sizeof(*new->S));
	new->U = calloc(m * n, sizeof(*new->U));
	new->VT = calloc(n * n, sizeof(*new->VT));
	complex double *work = calloc(1, sizeof(*work));
	double *rwork = calloc(lrwork, sizeof(*rwork));
	int *iwork = calloc(8 * n, sizeof(*iwork));

	if(!new->A || !new->A_new || !new->S || !new->U || !new->VT || !work || !rwork || !iwork) {
		SVD_free(new);
		new = NULL;
		goto done;
	}

	memcpy(new->A, A, m * n * sizeof(*new->A));
	memcpy(new->A_new, A, m * n * sizeof(*new->A_new));
	transpose(m, n, &new->A_new);

	LAPACK_zgesdd(&job, &m, &n, new->A_new, &lda, new->S, new->U, &ldu, new->VT, &ldvt, work, &lwork, rwork, iwork, &info);

	lwork = work[0];
	free(work);
	work = calloc(lwork, sizeof(*work));
	if(!work) {
		SVD_free(new);
		new = NULL;
		goto done;
	}

	LAPACK_zgesdd(&job, &m, &n, new->A_new, &lda, new->S, new->U, &ldu, new->VT, &ldvt, work, &lwork, rwork, iwork, &info);

	transpose(n, m, &new->A_new);
	transpose(n, m, &new->U);
	transpose(n, n, &new->VT);

	new->m = m;
	new->n = n;

done:
	free(work);
	free(rwork);
	free(iwork);
	return new;
}


#if 0
void print_matrix(char *filename, int rows, int cols, complex double *matrix)
{
	int i, j;
 	FILE *f = fopen(filename, "w");

 	for(i = 0; i < rows; i++) {
		for(j = 0; j < cols; j++)
			fprintf(f, "%g+I%g\t", creal(matrix[i * cols + j]), cimag(matrix[i * cols + j]));
		fprintf(f, "\n");
 	}
 	fclose(f);
}


void print_SVD_to_file(const struct SVD const *svd)
{
	int i;
	char filename[50];

	sprintf(filename, "radiometer_output/%s/S.dat", svd->SVD_name);
	FILE *S_file = fopen(filename, "w");
	for(i = 0; i < svd->n; i++)
		fprintf(S_file, "%g\n", svd->S[i]);
	fclose(S_file);

	sprintf(filename, "radiometer_output/%s/A.dat", svd->SVD_name);
	print_matrix(filename, svd->m, svd->n, svd->A);
	sprintf(filename, "radiometer_output/%s/A2.dat", svd->SVD_name);
	print_matrix(filename, svd->m, svd->n, svd->A_new);
	sprintf(filename, "radiometer_output/%s/U.dat", svd->SVD_name);
	print_matrix(filename, svd->m, svd->n, svd->U);
	sprintf(filename, "radiometer_output/%s/VT.dat", svd->SVD_name);
	print_matrix(filename, svd->n, svd->n, svd->VT);

	return;
}


void SVD_reverse_check(const struct SVD *svd)
{
	int i, j, k;
	char filename[50];
	double orig_row_mag;
	double reconst_row_mag;
	double diff_row_mag;

	complex double *US = calloc(svd->m * svd->n, sizeof(*US));
	complex double *USVT = calloc(svd->m * svd->n, sizeof(*USVT));

	if(!US || !USVT) {
		free(US);
		free(USVT);
		return;
	}

	for(i = 0; i < svd->m; i++)
		for(j = 0; j < svd->n; j++) {
			complex double value = 0;
			for(k = 0; k < svd->n; k++)
				value += US[i * svd->n + k] * svd->S[k] * svd->VT[k * svd->n + j];
			USVT[i * svd->n + j] = value;
			US[i * svd->n + j] = svd->U[i * svd->n + j] * svd->S[j];
		}

	sprintf(filename, "radiometer_output/%s/US.dat", svd->SVD_name);
	print_matrix(filename, svd->m, svd->n, US);
	sprintf(filename, "radiometer_output/%s/USVT.dat", svd->SVD_name);
	print_matrix(filename, svd->m, svd->n, USVT);


	sprintf(filename, "radiometer_output/%s/check1.dat", svd->SVD_name);
	FILE *error1 = fopen(filename, "w");
	sprintf(filename, "radiometer_output/%s/check2.dat", svd->SVD_name);
	FILE *error2 = fopen(filename, "w");

	for(i = 0; i < 1; i++) {
//	for(i = 0; i < svd->m; i++) {
//		reconst_row_mag = 0;
//		orig_row_mag = 0;
//		diff_row_mag = 0;
		for(j = 0; j < svd->n; j++) {
			complex double value = 0;
			for(k = 0; k < svd->n; k++) {
				fprintf(error1, "(%g %gI)\t*\t%g\t*\t(%g %gI)\t=\t%g %gI\n", creal(svd->U[i * svd->n + k]), cimag(svd->U[i * svd->n + k]), svd->S[k], creal(svd->VT[k * svd->n + j]), cimag(svd->VT[k * svd->n + j]), creal(svd->U[i * svd->n + k] * svd->S[k] * svd->VT[k * svd->n + j]), cimag(svd->U[i * svd->n + k] * svd->S[k] * svd->VT[k * svd->n + j]));
				value += svd->U[i * svd->n + k] * svd->S[k] * svd->VT[k * svd->n + j];
//				fprintf(error2, "value: %g %gI\n\n", creal(value), cimag(value));
			}
//			reconst_row_mag += value * conj(value);
//			orig_row_mag += svd->A[i * svd->n + j] * conj(svd->A[i * svd->n + j]);
//			diff_row_mag += (value - svd->A[i * svd->n + j]) * conj(value - svd->A[i * svd->n + j]);
//			fprintf(error2, "\n");
			fprintf(error1, "%g %g\n", creal(value), cimag(value));
		}
//		fprintf(error2, "%g\t%g\t%g\t%g\n", reconst_row_mag, orig_row_mag, diff_row_mag, diff_row_mag / orig_row_mag);
		fprintf(error1, "\n");
	}

	free(US);
	free(USVT);
	fclose(error1);
	fclose(error2);
	return;
}



// create plot data of VT
int SVD_to_sh_series_file(const struct SVD *svd, int l_max, int polar, double threshold)
{
	int i, j;
	struct sh_series *VT_sh = sh_series_new(l_max, polar);
	FILE *f = NULL;
	char filename[50];

	//for(i = 0; i < 20; i++) {
	for(i = 0; i < svd->n && svd->S[i] > threshold; i++) {
		for(j = 0; j < svd->n; j++) {
			VT_sh->coeff[j] = svd->VT[i * svd->n + j];
		}
		sprintf(filename, "radiometer_output/%s/dat_files/mode_%i.dat", svd->SVD_name, i);
		f = fopen(filename, "w");
		sh_series_print(f, VT_sh);
		fclose(f);
	}
	sh_series_free(VT_sh);
	return i;
}
#endif


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
	/* coefficients are in the same order in the sh_series_array as
	 * needed for A */
	struct SVD *first_SVD = SVD_new(array->n, sh_series_length(array->l_max, array->polar), array->coeff);

	return NULL;	/* FIXME */
}


#if 0
void new_sing(const struct sh_series_array *matrix)
{
	int i, j, l, k;

	/* coefficients are in the same order in the sh_series_array as
	 * needed for A */
	struct SVD *first_SVD = SVD_new(m, n, matrix->coeff);
	print_SVD_to_file(first_SVD);
	SVD_reverse_check(first_SVD);

	// construct new matrix to SVD
	int total_modes = SVD_to_sh_series_file(first_SVD, matrix->l_max, matrix->polar, 1);
	int total_phases = 2 * total_modes;
	complex double *A = calloc(total_phases * total_modes * n, sizeof(*A));
	if(!A) {
		free(A);
		return;
	}

	// due to multiple uses of m, k will be the m in (l,m) pairs
	for(i = 0; i < total_phases; i++)
		for(j = 0; j < total_modes; j++)
			for(l = 0; l <= l_max; l++)
				for(k = -l; k <= l; k++)
					A[i * total_modes * n + j * n + sh_series_params_lmoffset(l_max, l, k)] = cexp(2 * M_PI * I * i * k / total_phases) * first_SVD->VT[j * n + sh_series_params_lmoffset(l_max, l, k)] * first_SVD->S[j];

	struct SVD *second_SVD = SVD_new(total_phases * total_modes, n, A);
	print_SVD_to_file(second_SVD);
	SVD_reverse_check(second_SVD);
	SVD_to_sh_series_file(second_SVD, matrix->l_max, matrix->polar, 0);

	SVD_free(first_SVD);
	SVD_free(second_SVD);
	return;
}
#endif
