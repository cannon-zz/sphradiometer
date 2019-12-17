/*
 * Copyright (C) 2006--2009,2012,2019  Kipp C. Cannon
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


#ifndef __RADIOMETER_BIN_DIAGNOSTICS_H__
#define __RADIOMETER_BIN_DIAGNOSTICS_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <stdio.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/correlator.h>


#ifdef __cplusplus
extern "C" {
#define complex _Complex
#endif


/*
 * ============================================================================
 *
 *                                 Prototypes
 *
 * ============================================================================
 */


double randrange(double, double);
struct sh_series *random_sh_series(int, int);
void diagnostics_dump_sh_series(const struct sh_series *, char *);
int sh_series_cmp(const struct sh_series *, const struct sh_series *);
double diagnostics_rms_error(const struct sh_series *, const struct sh_series *);
void diagnostics_dump_correlator_plan_td_stats(FILE *, const struct correlator_plan_td *);
void diagnostics_dump_correlator_plan_fd_stats(FILE *, const struct correlator_plan_fd *);
void diagnostics_dump_network_plan_td_stats(FILE *, const struct correlator_network_plan_td *);
void diagnostics_dump_network_plan_fd_stats(FILE *, const struct correlator_network_plan_fd *);


#ifdef __cplusplus
#undef complex
}
#endif


#endif /* __RADIOMETER_BIN_DIAGNOSTICS_H__ */
