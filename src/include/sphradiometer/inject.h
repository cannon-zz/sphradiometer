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


#ifndef __RADIOMETER_INJECT_H__
#define __RADIOMETER_INJECT_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <gsl/gsl_rng.h>
#include <sphradiometer/instrument.h>


/*
 * ============================================================================
 *
 *                                 Prototypes
 *
 * ============================================================================
 */


double *gaussian_white_noise(double *, int, double, gsl_rng *);
double *sinusoidal_noise(double *, int, double, double, double, double);


double *inject_transform_point_source(const struct instrument *, double *, const double *, int, double, double, double);


double *inject_into(const struct instrument *, double *, const double *, int, double, double, double, double);


#endif  /* __RADIOMETER_INJECT_H__ */
