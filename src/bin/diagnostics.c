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


#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/correlator.h>
#include <diagnostics.h>


/*
 * ============================================================================
 *
 *                                Diagnostics
 *
 * ============================================================================
 */


/*
 * Convenience wrapper for printing an sh_series' coefficients to a file.
 */


void diagnostics_dump_sh_series(const struct sh_series *series, char *name)
{
	FILE *f = fopen(name, "w");

	sh_series_print(f, series);

	fclose(f);
}


/*
 * Compute the RMS difference between two functions on the sphere.
 */


double diagnostics_rms_error(const struct sh_series *measured, const struct sh_series *exact)
{
	double rms_error;

	/* determine the order of the difference */
	unsigned int l = (measured->l_max > exact->l_max) ? measured->l_max : exact->l_max;

	/* make working copies */
	struct sh_series *m_copy = sh_series_copy(measured);
	struct sh_series *e_copy = sh_series_copy(exact);

	/* bring the working copies to a common size */
	m_copy = sh_series_resize(m_copy, l);
	e_copy = sh_series_resize(e_copy, l);
	if(!(measured->polar && exact->polar)) {
		sh_series_set_polar(m_copy, 0);
		sh_series_set_polar(e_copy, 0);
	}

	/* compute the difference */
	sh_series_sub(e_copy, 1.0, m_copy);

	/* compute the RMS of the difference */
	rms_error = sqrt(creal(sh_series_dot(e_copy, e_copy)) / (4 * M_PI));

	/* clean up */
	sh_series_free(m_copy);
	sh_series_free(e_copy);

	return rms_error;
}


/*
 * Print correlator and baseline summary information.
 */


void diagnostics_dump_correlator_plan_td_stats(FILE *f, const struct correlator_plan_td *plan)
{
	sh_series_print(f, &plan->proj_a->series[(plan->proj_a->n - 1) / 2]);
	{
	gsl_vector *phase_centre = instrument_array_get(plan->baseline->instruments, plan->baseline->index_a)->phase_centre;
	fprintf(f, "instrument A: phase centre = (%g, %g, %g)\n", gsl_vector_get(phase_centre, 0), gsl_vector_get(phase_centre, 1), gsl_vector_get(phase_centre, 2));
	}
	fprintf(f, "instrument A: n = %d\n", plan->proj_a->n);
	fprintf(f, "instrument A: l max = %d\n", plan->proj_a->l_max);
	{
	gsl_vector *phase_centre = instrument_array_get(plan->baseline->instruments, plan->baseline->index_b)->phase_centre;
	fprintf(f, "instrument B: phase centre = (%g, %g, %g)\n", gsl_vector_get(phase_centre, 0), gsl_vector_get(phase_centre, 1), gsl_vector_get(phase_centre, 2));
	}
	fprintf(f, "instrument B: n = %d\n", plan->proj_b->n);
	fprintf(f, "instrument B: l max = %d\n", plan->proj_b->l_max);
	fprintf(f, "correlator: transient = %d samples\n", plan->transient);
	fprintf(f, "correlator: l max = %d\n", plan->power_2d->l_max);
	fprintf(f, "correlator: azimuthal symmetry = %s\n", plan->power_1d->polar ? "enabled" : "disabled");
}


void diagnostics_dump_correlator_plan_fd_stats(FILE *f, const struct correlator_plan_fd *plan)
{
	{
	gsl_vector *phase_centre = instrument_array_get(plan->baseline->instruments, plan->baseline->index_a)->phase_centre;
	fprintf(f, "instrument A: phase centre = (%g, %g, %g)\n", gsl_vector_get(phase_centre, 0), gsl_vector_get(phase_centre, 1), gsl_vector_get(phase_centre, 2));
	}
	{
	gsl_vector *phase_centre = instrument_array_get(plan->baseline->instruments, plan->baseline->index_b)->phase_centre;
	fprintf(f, "instrument B: phase centre = (%g, %g, %g)\n", gsl_vector_get(phase_centre, 0), gsl_vector_get(phase_centre, 1), gsl_vector_get(phase_centre, 2));
	}
	fprintf(f, "correlator: N = %d samples\n", plan->delay_product->n);
	fprintf(f, "correlator: transient = %d samples\n", plan->transient);
	fprintf(f, "correlator: l max = %d\n", plan->power_2d->l_max);
	fprintf(f, "correlator: azimuthal symmetry = %s\n", plan->power_1d->polar ? "enabled" : "disabled");
}


void diagnostics_dump_network_plan_td_stats(FILE *f, const struct correlator_network_plan_td *plan)
{
	int i;
	unsigned l_max = 0;

	for(i = 0; i < plan->baselines->n_baselines; i++) {
		fprintf(f, "=== Baseline %d ===\n", i);
		diagnostics_dump_correlator_plan_td_stats(f, plan->plans[i]);
		if(plan->plans[i]->power_2d->l_max > l_max)
			l_max = plan->plans[i]->power_2d->l_max;
	}
	fprintf(f, "=== End Baselines ===\n");
	fprintf(f, "sky: l max = %d\n", l_max);
}


void diagnostics_dump_network_plan_fd_stats(FILE *f, const struct correlator_network_plan_fd *plan)
{
	int i;
	unsigned l_max = 0;

	for(i = 0; i < plan->baselines->n_baselines; i++) {
		fprintf(f, "=== Baseline %d ===\n", i);
		diagnostics_dump_correlator_plan_fd_stats(f, plan->plans[i]);
		if(plan->plans[i]->power_2d->l_max > l_max)
			l_max = plan->plans[i]->power_2d->l_max;
	}
	fprintf(f, "=== End Baselines ===\n");
	fprintf(f, "sky: l max = %d\n", l_max);
}
