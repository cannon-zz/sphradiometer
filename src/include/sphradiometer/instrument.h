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
#ifndef SWIG
#include <lal/LALDetectors.h>
#endif


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
	/*
	 * user data.  not used by the sphradiometer library.  if freefunc
	 * is not NULL it will be called with data as the one argument when
	 * the instrument object is free()ed, otherwise the memory at the
	 * location pointed to by data will not be free()ed when the
	 * instrument object is free()ed.
	 */
	void *data;
	void (*freefunc)(void *);
};


/*
 * An array of instruments
 */


struct instrument_array {
	struct instrument **instruments;
	int n;
};


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


struct instrument *instrument_new(double, double, double, void *, void (*)(void *));
void instrument_free(struct instrument *);
struct instrument *instrument_new_from_r_theta_phi(double, double, double, void *, void (*)(void *));
#ifndef SWIG
struct instrument *instrument_new_from_LALDetector(const LALDetector *);
#endif
struct instrument *instrument_new_from_name(const char *);


struct instrument_array *instrument_array_new(int);
int instrument_array_len(const struct instrument_array *);
struct instrument *instrument_array_get(const struct instrument_array *, int);
struct instrument *instrument_array_set(struct instrument_array *, int, struct instrument *);
struct instrument *instrument_array_append(struct instrument_array *, struct instrument *);
void instrument_array_free(struct instrument_array *);


double instrument_r(const struct instrument *);
gsl_vector *instrument_baseline(const struct instrument *, const struct instrument *);


#endif  /* __RADIOMETER_INSTRUMENT_H__ */
