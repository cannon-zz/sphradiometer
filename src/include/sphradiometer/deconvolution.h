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


#ifndef __RADIOMETER_DECONVOLUTION_H__
#define __RADIOMETER_DECONVOLUTION_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <sphradiometer/sh_series.h>


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


double complex *sh_series_array_orthogonalize(struct sh_series_array *);


#ifdef __cplusplus
#undef complex
}
#endif


#endif  /* __RADIOMETER_DECONVOLUTION_H__ */
