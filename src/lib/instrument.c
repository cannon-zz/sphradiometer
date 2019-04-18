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


#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <sphradiometer/instrument.h>


/*
 * ============================================================================
 *
 *                                BLAS Extras
 *
 * ============================================================================
 */


/*
 * BLAS-like routines for common tasks.
 */


/*
 * Compute the length of a 3-vector.
 */


double vector_magnitude(const gsl_vector *v)
{
	const double x = gsl_vector_get(v, 0);
	const double y = gsl_vector_get(v, 1);
	const double z = gsl_vector_get(v, 2);

	return sqrt(x * x + y * y + z * z);
}


/*
 * Compute inner product of 3-vector \vec{v} with the unit vector \hat{s}
 * pointing in the direction given by the spherical polar co-ordinates
 * theta, phi.
 */


double vector_r_dot_s(const gsl_vector *v, double theta, double phi)
{
	return (gsl_vector_get(v, 0) * cos(phi) + gsl_vector_get(v, 1) * sin(phi)) * sin(theta) + gsl_vector_get(v, 2) * cos(theta);
}


/*
 * Compute the spherical polar co-ordinates, theta and phi, describing the
 * direction in which a 3-vector points.  theta is in [0, pi) and phi will
 * be in [0, 2 pi)
 */


void vector_direction(const gsl_vector *v, double *theta, double *phi)
{
	const double x = gsl_vector_get(v, 0);
	const double y = gsl_vector_get(v, 1);
	const double z = gsl_vector_get(v, 2);

	*theta = atan2(sqrt(x * x + y * y), z);
	*phi = atan2(y, x);
	if(*phi < 0)
		*phi += 2 * M_PI;
}


/*
 * ============================================================================
 *
 *                             Instrument Object
 *
 * ============================================================================
 */


/*
 * Create and destroy an instrument object.
 */


struct instrument *instrument_new(double x, double y, double z, void *data, void (*freefunc)(void *))
{
	struct instrument *new = malloc(sizeof(*new));
	gsl_vector *phase_centre = gsl_vector_alloc(3);

	if(!new || !phase_centre) {
		free(new);
		gsl_vector_free(phase_centre);
		if(freefunc)
			freefunc(data);
		return NULL;
	}

	gsl_vector_set(phase_centre, 0, x);
	gsl_vector_set(phase_centre, 1, y);
	gsl_vector_set(phase_centre, 2, z);

	new->phase_centre = phase_centre;
	new->data = data;
	new->freefunc = freefunc;

	return new;
}


void instrument_free(struct instrument *instrument)
{
	if(instrument) {
		gsl_vector_free(instrument->phase_centre);
		if(instrument->freefunc)
			instrument->freefunc(instrument->data);
	}
	free(instrument);
}


struct instrument *instrument_new_from_r_theta_phi(double r, double theta, double phi, void *data, void (*freefunc)(void *))
{
	return instrument_new(
		r * sin(theta) * cos(phi),
		r * sin(theta) * sin(phi),
		r * cos(theta),
		data,
		freefunc
	);
}


/*
 * Distance from origin to instrument phase centre.
 */


double instrument_r(const struct instrument *instrument)
{
	return vector_magnitude(instrument->phase_centre);
}


/*
 * Vector separation of two instruments' phase centres.
 * instrument_baseline(a,b) returns the vector pointing from b's phase
 * centre to a's:  like a subtraction operator, returns a - b.
 */


gsl_vector *instrument_baseline(const struct instrument *a, const struct instrument *b)
{
	gsl_vector *baseline = gsl_vector_alloc(3);

	gsl_vector_memcpy(baseline, a->phase_centre);
	gsl_vector_sub(baseline, b->phase_centre);

	return baseline;
}


/*
 * ============================================================================
 *
 *                          Instrument Array Object
 *
 * ============================================================================
 */


struct instrument_array *instrument_array_new(int n)
{
	struct instrument_array *instrument_array = malloc(sizeof(*instrument_array));
	struct instrument **instruments = calloc(n, sizeof(*instruments));

	if(!instrument_array || !instruments) {
		free(instrument_array);
		free(instruments);
		return NULL;
	}

	instrument_array->instruments = instruments;
	instrument_array->n = n;

	return instrument_array;
}


int instrument_array_len(const struct instrument_array *instrument_array)
{
	return instrument_array->n;
}


/*
 * returns borrowed reference
 */

struct instrument *instrument_array_get(const struct instrument_array *instrument_array, int k)
{
	return 0 <= k && k < instrument_array->n ? instrument_array->instruments[k] : NULL;
}


/*
 * takes ownership of instrument;  frees it on failure
 */

struct instrument *instrument_array_set(struct instrument_array *instrument_array, int k, struct instrument *instrument)
{
	if(k < 0 || k >= instrument_array->n) {
		instrument_free(instrument);
		NULL;
	}

	instrument_free(instrument_array->instruments[k]);
	instrument_array->instruments[k] = instrument;
	return instrument;
}


/*
 * takes ownership of instrument;  frees it on failure
 */


struct instrument *instrument_array_append(struct instrument_array *instrument_array, struct instrument *instrument)
{
	struct instrument **new = realloc(instrument_array->instruments, (instrument_array->n + 1) * sizeof(*instrument_array->instruments));

	if(!new) {
		instrument_free(instrument);
		return NULL;
	}

	instrument_array->instruments = new;
	instrument_array->instruments[instrument_array->n] = instrument;
	instrument_array->n++;

	return instrument;
}


/*
 * frees instrument_array and all instruments it contains
 */


void instrument_array_free(struct instrument_array *instrument_array)
{
	if(instrument_array) {
		int i;
		for(i = 0; i < instrument_array->n; i++)
			instrument_free(instrument_array->instruments[i]);
		free(instrument_array->instruments);
	}
	free(instrument_array);
}
