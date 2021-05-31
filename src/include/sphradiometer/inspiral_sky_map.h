/*
 * Copyright (C) 2020  Takuya Tsutsui
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


#include <sphradiometer/correlator.h>
#include <sphradiometer/sh_series.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>


/*
 * ============================================================================
 *
 *                                    API
 *
 * ============================================================================
 */


struct sh_series *sh_series_log_uniformsky_prior(int);

int correlator_network_plan_mult_by_projection(struct correlator_network_plan_fd *, double, double, double **);

int generate_alm_skys(struct sh_series **, struct sh_series **, struct correlator_network_plan_fd *, struct correlator_network_plan_fd *, COMPLEX16TimeSeries **, COMPLEX16Sequence **, struct sh_series *);
