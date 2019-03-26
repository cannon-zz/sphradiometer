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


#ifndef __RADIOMETER_BIN_INSTRUMENTS_H__
#define __RADIOMETER_BIN_INSTRUMENTS_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <sphradiometer/instrument.h>
#include <lal/LALDetectors.h>


/*
 * ============================================================================
 *
 *                                 Prototypes
 *
 * ============================================================================
 */


struct instrument *instrument_from_r_theta_phi(double, double, double, void *, void (*)(void *));
struct instrument *instrument_from_LALDetector(const LALDetector *);
struct instrument *instrument_from_name(const char *);


#endif /* __RADIOMETER_BIN_INSTRUMENTS_H__ */
