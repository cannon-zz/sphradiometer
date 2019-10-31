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


#include <complex.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <sphradiometer/instrument.h>


/*
 * ============================================================================
 *
 *                              Noise Generators
 *
 * ============================================================================
 */


double *gaussian_white_noise(double *series, int n, double variance, gsl_rng *rng)
{
	int i;

	for(i = 0; i < n; i++)
		series[i] = gsl_ran_gaussian(rng, sqrt(variance));

	return series;
}


double *sinusoidal_noise(double *series, int n, double delta_t, double variance, double f, double phase)
{
	const double amplitude = sqrt(2 * variance);
	const double omega = 2 * M_PI * delta_t * f;
	int i;

	for(i = 0; i < n; i++)
		series[i] = amplitude * sin(omega * i + phase);

	return series;
}


/*
 * ============================================================================
 *
 *                     Sky --> Instrument Transformation
 *
 * ============================================================================
 */


/*
 * Transform a time series for injection into the given instrument's
 * output, assuming the time series originates from a point source at the
 * given location.
 */


static double *_inject_transform_point_source(const struct instrument *instrument, double *output, const double *input, const int n, const double delta_t, const double theta, const double phi)
{
	/* This implementation does not work for noise, only for nice
	 * smooth functions that are well over-sampled.  I'm leaving it
	 * here for reference purposes. */
	double t[n];
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_cspline, n);
	gsl_interp_accel *accel = gsl_interp_accel_alloc();
	const double r_dot_s = vector_r_dot_s(instrument->phase_centre, theta, phi);
	const int transient = ceil(fabs(r_dot_s / delta_t)) + 1;
	int i;

	/* generate values for the time co-ordinate */
	for(i = 0; i < n; i++)
		t[i] = i * delta_t;

	memset(output, 0, transient * sizeof(*output));
	memset(output + n - transient, 0, transient * sizeof(*output));
	for(i = r_dot_s < 0.0 ? transient : 0; i < n - (r_dot_s > 0.0 ? transient : 0); i++)
		output[i] = gsl_interp_eval(interp, t, input, i * delta_t + r_dot_s, accel);

	gsl_interp_accel_free(accel);
	gsl_interp_free(interp);

	return output;
}


double *inject_transform_point_source(const struct instrument *instrument, double *output, const double *input, const int n, double const delta_t, const double theta, const double phi)
{
	/* This works very well even for noise, but treats the data as
	 * periodic, wrapping it around the end.  If this is inappropriate,
	 * the data should be windowed in advance. */
	complex double *ft = malloc((n / 2 + 1) * sizeof(*ft));
	fftw_plan fwd, rev;
	const double r_dot_s = vector_r_dot_s(instrument->phase_centre, theta, phi);
	const complex double w = I * 2 * M_PI * r_dot_s / (delta_t * n);
	int i;

	if(!ft)
		return NULL;

	fwd = fftw_plan_dft_r2c_1d(n, input, ft, FFTW_ESTIMATE);
	rev = fftw_plan_dft_c2r_1d(n, ft, output, FFTW_ESTIMATE);
	if(!fwd || !rev) {
		fftw_destroy_plan(rev);
		fftw_destroy_plan(fwd);
		free(ft);
		return NULL;
	}

	fftw_execute(fwd);
	for(i = 0; i < n / 2 + 1; i++)
		ft[i] *= cexp(w * i) / n;
	fftw_execute(rev);

	fftw_destroy_plan(rev);
	fftw_destroy_plan(fwd);
	free(ft);

	return output;
}


/*
 * ============================================================================
 *
 *                                  Wrapper
 *
 * ============================================================================
 */


/*
 * Accept an instrument, a time series that is the instrument's nominal
 * output, a time series for an injection to add to the instrument's
 * output, the length of the series and the sample rate, the right
 * ascension and declination of the source of the injection, and the GMST,
 * and transform the injection time series appropriately and add it to the
 * instrument's output.
 *
 * Note that, for now, this function assumes the instrument is stationary
 * for the duration of the injection, and so can only be used to inject
 * signals for which this assumption is appropriate (inspirals, ring-downs,
 * super nova bursts, etc.).  Don't use this to inject hours and hours of a
 * pulsar signal.
 *
 * Also note that this function does not know anything about antenna
 * patterns or frequency responses.  Just geometric delays.
 *
 * Basically it's pretty useless except for simple demos.
 */


double *inject_into(const struct instrument *instrument, double *series, const double *injection, int n, double delta_t, double ra, double dec, double gmst)
{
	double *transformed_injection = malloc(n * sizeof(*transformed_injection));
	const double greenwich_hour_angle = gmst - ra;
	int i;

	if(!transformed_injection)
		return NULL;

	inject_transform_point_source(instrument, transformed_injection, injection, n, delta_t, M_PI_2 - dec, -greenwich_hour_angle);
	for(i = 0; i < n; i++)
		series[i] += transformed_injection[i];

	free(transformed_injection);

	return series;
}
