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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <sphradiometer/sky.h>


/*
 * ============================================================================
 *
 *                                 Rotations
 *
 * ============================================================================
 */


/*
 * The matrix returned by sky_rotation_matrix() rotates by alpha about the
 * positive z axis, then by beta about the positive y axis, then by gamma
 * about the positive z axis.  sky_inv_rotation_matrix() returns the
 * inverse of the same matrix.
 */


double *sky_rotation_matrix(double alpha, double beta, double gamma)
{
	enum {
		x = 0,
		y = 1,
		z = 2
	};
	double *R = malloc(9 * sizeof(*R));

	if(!R)
		return NULL;

	R[3 * x + x] =  cos(gamma) * cos(beta) * cos(alpha) - sin(gamma) * sin(alpha);
	R[3 * x + y] = -cos(gamma) * cos(beta) * sin(alpha) - sin(gamma) * cos(alpha);
	R[3 * x + z] =  cos(gamma) * sin(beta);
	R[3 * y + x] =  sin(gamma) * cos(beta) * cos(alpha) - cos(gamma) * sin(alpha);
	R[3 * y + y] = -sin(gamma) * cos(beta) * sin(alpha) + cos(gamma) * cos(alpha);
	R[3 * y + z] =  sin(gamma) * sin(beta);
	R[3 * z + x] = -sin(beta) * cos(alpha);
	R[3 * z + y] =  sin(beta) * sin(alpha);
	R[3 * z + z] =  cos(beta);

	return R;
}


double *sky_inv_rotation_matrix(double alpha, double beta, double gamma)
{
	return sky_rotation_matrix(-gamma, -beta, -alpha);
}


/*
 * The matrix returned by sh_series_equatorial_to_galactic_rot_matrix()
 * rotates the right-handed equatorial co-ordinate system whose x axis
 * points to the vernal equinox and whose z axis points to the Earth's
 * north celestial pole to the right-handed galactic co-ordinate system
 * whose x axis points to the galactic core and whose z axis points to the
 * galactic pole.  sky_galactic_to_equatorial_rot_matrix() returns the
 * inverse of this matrix.
 */


double *sky_equatorial_to_galactic_rot_matrix(void)
{
	return sky_rotation_matrix(-SKY_MW_Z_RA_RAD, M_PI_2 - SKY_MW_Z_DEC_RAD, SKY_MW_Z_RA_RAD - acos(cos(SKY_MW_X_J2000_DEC_RAD) * cos(SKY_MW_X_J2000_RA_RAD)));
}


double *sky_galactic_to_equatorial_rot_matrix(void)
{
	return sky_rotation_matrix(SKY_MW_Z_RA_RAD - acos(cos(SKY_MW_X_J2000_DEC_RAD) * cos(SKY_MW_X_J2000_RA_RAD)), M_PI_2 - SKY_MW_Z_DEC_RAD, -SKY_MW_Z_RA_RAD);
}
