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


#ifndef __RADIOMETER_INSTRUMENT_H__
#define __RADIOMETER_INSTRUMENT_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <gsl/gsl_vector.h>


/*
 * ============================================================================
 *
 *                                 Data Types
 *
 * ============================================================================
 */


/*
 * Description of one element in the detector network.
 */


struct instrument {
	gsl_vector *phase_centre;
};


/*
 * ============================================================================
 *
 *                                   Macros
 *
 * ============================================================================
 */


/*
 * ============================================================================
 *
 *                                 Prototypes
 *
 * ============================================================================
 */


double vector_magnitude(const gsl_vector *);
double vector_r_dot_s(const gsl_vector *, double, double);
void vector_direction(const gsl_vector *, double *, double *);


struct instrument *instrument_new(double, double, double);
void instrument_free(struct instrument *);


double instrument_r(const struct instrument *);
gsl_vector *instrument_baseline(const struct instrument *, const struct instrument *);


#endif  /* __RADIOMETER_INSTRUMENT_H__ */
