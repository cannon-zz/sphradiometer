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


#include <stdlib.h>
#include <gsl/gsl_const.h>
#include <sphradiometer/instrument.h>
#include <lal/LALDetectors.h>
#include <lal/LALSimulation.h>
#include <instruments.h>


/*
 * ============================================================================
 *
 *                                Instruments
 *
 * ============================================================================
 */


struct instrument *instrument_from_LALDetector(const LALDetector *det)
{
	/* this copy is only being made to remove the const'edness */
	LALDetector *copy = malloc(sizeof(*copy));
	if(!copy)
		return NULL;
	*copy = *det;
	return instrument_new(
		det->location[0] / GSL_CONST_MKS_SPEED_OF_LIGHT,
		det->location[1] / GSL_CONST_MKS_SPEED_OF_LIGHT,
		det->location[2] / GSL_CONST_MKS_SPEED_OF_LIGHT,
		copy,
		free
	);
}


struct instrument *instrument_from_name(const char *name)
{
	return instrument_from_LALDetector(XLALDetectorPrefixToLALDetector(name));
}
