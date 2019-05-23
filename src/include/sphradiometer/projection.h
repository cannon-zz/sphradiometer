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


#ifndef __RADIOMETER_PROJECTION_H__
#define __RADIOMETER_PROJECTION_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <gsl/gsl_vector.h>
#include <sphradiometer/sh_series.h>


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
 * Aggregate data passed to delay function while projecting onto spherical
 * harmonics.
 */


struct projection_delay_element_data {
	const gsl_vector *r;
	double delta_t;
	int j, k;
};


/*
 * ============================================================================
 *
 *                                 Prototypes
 *
 * ============================================================================
 */


int projection_matrix_n_elements(double , double);
unsigned int projection_matrix_l_max(double, double);


double projection_delay_element(double, double, void *);
struct sh_series_array *projection_matrix_delay(unsigned int, unsigned int, const gsl_vector *, double);


#ifdef __cplusplus
#undef complex
}
#endif


#endif  /* __RADIOMETER_PROJECTION_H__ */
