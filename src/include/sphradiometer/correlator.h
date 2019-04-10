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


#ifndef __RADIOMETER_CORRELATOR_H__
#define __RADIOMETER_CORRELATOR_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <fftw3.h>
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                 Data Types
 *
 * ============================================================================
 */


/*
 * Description of a correlator baseline:  a pair of instruments, with some
 * pre-computed metadata about their geometry.
 */


struct correlator_baseline {
	const struct instrument_array *instruments;
	int index_a;
	int index_b;
	gsl_vector *d;
	double theta;
	double phi;
};


/*
 * Time-domain correlation plan, a collection of pre-computed information
 * and pre-allocated storage for use in cross-correlating data for a
 * baseline.  This is specficially for the purely time-domain case, in
 * which the input vectors contain time series data.  This case allows the
 * vector sizes to go unspecified, and the correlator will accept any
 * amount of input data.
 */


struct correlator_plan_td {
	const struct correlator_baseline *baseline;
	double delta_t;
	int transient;
	struct sh_series_array *proj_a;
	struct sh_series_array *proj_b;
	struct sh_series *sample_a;
	struct sh_series *sample_b;
	struct sh_series *product;
	struct sh_series *power_1d;
	struct sh_series *power_2d;
	struct sh_series_product_plan *product_plan;
	struct sh_series_rotation_plan *rotation_plan;
};


/*
 * Frequency-domain correlation plan.  If the input vector sizes are known
 * in advance, the frequency-domain correlation technique can be used.  The
 * frequency-domain correlator is much faster than the time-domain
 * correlator.  This is due to the ability to pre-compute the spherical
 * harmonic series products.
 *
 * It is simply the advance knowledge of the input vector's size that
 * allows the SH products to be pre-computed, so this could be done in the
 * time-domain case as well.  However, if one is willing to sacrifice the
 * freedom of having arbitrary input vector sizes, then one might as well
 * Fourier transform the input vectors and also get the speed advantage of
 * frequency-domain convolution.  In this way, there is no reason to
 * implement the pre-computed products in the time-domain case.
 */


struct correlator_plan_fd {
	const struct correlator_baseline *baseline;
	double delta_t;
	int transient;
	struct sh_series_array *delay_product;
	complex double *fseries_product;
	struct sh_series *power_1d;
	struct sh_series *power_2d;
	struct sh_series_rotation_plan *rotation_plan;
};


/*
 * A network of baselines.
 */


struct correlator_network_baselines {
	struct correlator_baseline **baselines;
	int n_baselines;
};


/*
 * Time-domain correlation plan for a baseline network.
 */


struct correlator_network_plan_td {
	struct correlator_network_baselines *baselines;
	struct correlator_plan_td **plans;
};


/*
 * Frequency-domain correlation plan for a baseline network.
 */


struct correlator_network_plan_fd {
	struct correlator_network_baselines *baselines;
	struct correlator_plan_fd **plans;
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


/*
 * Utilties
 */


unsigned int correlator_baseline_power_l_max(const struct correlator_baseline *, double);
double correlator_dump_interval(unsigned int, unsigned int);
int correlator_transient(const struct sh_series_array *, const struct sh_series_array *);


/*
 * Windowing
 */


double *correlator_tukey_window_new(int, int, int, double);
double *correlator_square_window_new(int, int, double);


/*
 * Fourier transforms
 */


fftw_plan correlator_tseries_to_fseries_plan(double *, complex double *, int);
void correlator_tseries_to_fseries(double *, complex double *, int, fftw_plan);

fftw_plan correlator_ctseries_to_fseries_plan(complex double *, complex double *, int);
void correlator_ctseries_to_fseries(fftw_plan);


/*
 * Baselines
 */


struct correlator_baseline *correlator_baseline_new(const struct instrument_array *, int, int);
void correlator_baseline_free(struct correlator_baseline *);


struct correlator_network_baselines *correlator_network_baselines_new(const struct instrument_array *);
void correlator_network_baselines_free(struct correlator_network_baselines *);
unsigned int correlator_network_l_max(struct correlator_network_baselines *, double);


/*
 * Correlation plans
 */


struct correlator_plan_td *correlator_plan_td_new(const struct correlator_baseline *, double);
void correlator_plan_td_free(struct correlator_plan_td *);


struct correlator_plan_fd *correlator_plan_fd_new(const struct correlator_baseline *, int, double);
void correlator_plan_fd_free(struct correlator_plan_fd *);


struct correlator_network_plan_td *correlator_network_plan_td_new(struct correlator_network_baselines *, double);
void correlator_network_plan_td_free(struct correlator_network_plan_td *);


struct correlator_network_plan_fd *correlator_network_plan_fd_new(struct correlator_network_baselines *, int, double);
void correlator_network_plan_fd_free(struct correlator_network_plan_fd *);


/*
 * Correlators
 */


struct sh_series *correlator_baseline_integrate_power_td(const double *, const double *, const double *, int, struct correlator_plan_td *);


struct sh_series *correlator_baseline_integrate_power_fd(const complex double *, const complex double *, struct correlator_plan_fd *);


struct sh_series *correlator_network_integrate_power_td(struct sh_series *, double **, int, double **, struct correlator_network_plan_td *);


struct sh_series *correlator_network_integrate_power_fd(struct sh_series *, complex double **, struct correlator_network_plan_fd *);


#endif  /* __RADIOMETER_CORRELATOR_H__ */
