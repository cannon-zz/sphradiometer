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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/projection.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                          projection_matrix Object
 *
 * ============================================================================
 */


/*
 * Desired number of elements in projection matrix based on the sample rate
 * and the distance from the origin to the phase centre.  The result is
 * guaranteed to be odd.
 */


int projection_matrix_n_elements(double r, double delta_t)
{
#if 1
	const int buffer = 60;
#else
	/* alternative for doing tests in Phys. Rev. paper */
	const int buffer = 0;
#endif

	return 2 * (int) (buffer + r / delta_t + 1) + 1;
}


/*
 * Desired spherical harmonic order for projection matrix elements based on
 * the sample rate and the distance from the origin to the phase centre.
 */


unsigned int projection_matrix_l_max(double r, double delta_t)
{
	const int l_max = M_PI * r / delta_t + 1;

#if 1
	return l_max + 4;
#else
	/* alternative for doing tests in Phys. Rev. paper */
	return l_max;
#endif
}


/*
 * ============================================================================
 *
 *                     Instrument --> Sky Transformation
 *
 * ============================================================================
 */


double projection_delay_element(double theta, double phi, void *data)
{
	const double r_dot_s = vector_r_dot_s(((struct projection_delay_element_data *) data)->r, theta, phi);
	const double delta_t = ((struct projection_delay_element_data *) data)->delta_t;
	const int j = ((struct projection_delay_element_data *) data)->j;
	const int k = ((struct projection_delay_element_data *) data)->k;
	const double pi_x = M_PI * (j - k - r_dot_s / delta_t);

	return (pi_x == 0.0) ? 1.0 : sin(pi_x) / pi_x;
}


struct sh_series_array *projection_matrix_delay(unsigned int n, unsigned int l_max, const gsl_vector *phase_centre, double delta_t)
{
	const int polar = (gsl_vector_get(phase_centre, 0) == 0.0) && (gsl_vector_get(phase_centre, 1) == 0.0);
	struct sh_series_array *matrix = sh_series_array_new(n, l_max, polar);
	struct projection_delay_element_data data = {
		.r = phase_centre,
		.delta_t = delta_t,
		.j = 0
	};
	int i;

	if(!matrix)
		return NULL;

	for(i = 0, data.k = -(int) (n - 1) / 2; i < (int) n; i++, data.k++)
		sh_series_from_realfunc(&matrix->series[i], projection_delay_element, &data);

	return matrix;
}


/*
 * ============================================================================
 *
 *                                  File I/O
 *
 * ============================================================================
 */


/*
 * Write a projection_matrix object to a file.  Warning:  the data is not
 * portable;  a projection_matrix object can only be reconstructed from the
 * file on a machine with the same word size, alignment, and endianness.
 * Returns 0 on success, -1 on failure (use ferror() or errno for details).
 *
 * This is a hack, don't use in production code.
 */


int projection_matrix_write(const struct sh_series_array *matrix, FILE *file)
{
	int i;
	int n;

	n = fwrite(&matrix->n, sizeof(matrix->n), 1, file);
	n += fwrite(&matrix->l_max, sizeof(matrix->l_max), 1, file);
	n += fwrite(&matrix->polar, sizeof(matrix->polar), 1, file);
	for(i = 0; i < matrix->n; i++)
		n += sh_series_write(&matrix->series[i], file) == 0;

	return -(n != 2 + matrix->n);
}


/*
 * Read a projection_matrix from a file.  Returns a pointer to the
 * projection matrix on success, NULL on failure (use ferror() or errno for
 * details).
 *
 * This is a hack, don't use in production code.
 */


struct sh_series_array *projection_matrix_read(FILE *file)
{
	struct sh_series_array *matrix;
	unsigned int n;
	unsigned int l_max;
	int polar;
	unsigned int i;

	fread(&n, sizeof(n), 1, file);
	fread(&l_max, sizeof(l_max), 1, file);
	fread(&polar, sizeof(polar), 1, file);
	matrix = sh_series_array_new(n, l_max, polar);
	for(i = 0; i < n; i++)
		sh_series_read(&matrix->series[i], file);

	return matrix;
}
