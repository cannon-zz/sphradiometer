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
#include <gsl/gsl_math.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/sky.h>


/*
 * ============================================================================
 *
 *                                 Rotations
 *
 * ============================================================================
 */


/*
 * The matrix returned by sky_equatorial_to_galactic_rot_matrix() rotates
 * the right-handed equatorial co-ordinate system whose x axis points to
 * the vernal equinox and whose z axis points to the Earth's north
 * celestial pole to the right-handed galactic co-ordinate system whose x
 * axis points to the galactic core and whose z axis points to the galactic
 * pole.  sky_galactic_to_equatorial_rot_matrix() returns the inverse of
 * this matrix.
 */


double *sky_equatorial_to_galactic_rot_matrix(void)
{
	/* the alpha parameter is the angle between the X axis and the line
	 * of nodes.  this is obtained knowing the co-ordinates of those
	 * two points and the law of spherical cosines
	 *
	 * the line of nodes is pi/2 from the positive Z axis, on a merdian
	 * at right ascension = (pi/2 + Z axis RA).  the X axis is (pi/2 -
	 * X axis DEC) from the positive Z axis, on a meridian at right
	 * ascension = (X axis RA).
	 */

	return euler_rotation_matrix(
		-SKY_MW_Z_J2000_RA_RAD,
		SKY_MW_Z_J2000_DEC_RAD - M_PI_2,
		asin(cos(SKY_MW_X_J2000_DEC_RAD) * sin(SKY_MW_Z_J2000_RA_RAD - SKY_MW_X_J2000_RA_RAD))
	);
}


double *sky_galactic_to_equatorial_rot_matrix(void)
{
	return euler_inv_rotation_matrix(
		-SKY_MW_Z_J2000_RA_RAD,
		SKY_MW_Z_J2000_DEC_RAD - M_PI_2,
		asin(cos(SKY_MW_X_J2000_DEC_RAD) * sin(SKY_MW_Z_J2000_RA_RAD - SKY_MW_X_J2000_RA_RAD))
	);
}
