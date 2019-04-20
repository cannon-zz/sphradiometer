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
#include <sphradiometer/sky.h>


/*
 * ============================================================================
 *
 *                                 Rotations
 *
 * ============================================================================
 */


/*
 * The matrix returned by euler_rotation_matrix() rotates by gamma about
 * the positive z axis, then by beta about the positive y axis, then by
 * alpha about the positive z axis.  euler_inv_rotation_matrix() returns
 * the inverse of the same matrix.
 *
 * Compared to wikipedia's article on Euler angles, this is the
 * "y-convention Geometrical definition", with our gamma = wikipedia's
 * alpha, our beta = wikipedia's beta, our alpha = wikipedia's gamma.
 *
 * https://en.wikipedia.org/wiki/Euler_angles
 */


double *euler_rotation_matrix(double gamma, double beta, double alpha)
{
	enum {
		x = 0,
		y = 1,
		z = 2
	};
	double *R = malloc(9 * sizeof(*R));

	if(!R)
		return NULL;

	/*
	 * using rotation matrixes from
	 *
	 * https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
	 *
	 * wxMaxima code to generate these expressions below:
	 *
	 * Rz(w) := matrix([cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]);
	 * Ry(w) := matrix([cos(w), 0, sin(w)], [0, 1, 0], [-sin(w), 0, cos(w)]);
	 *
	 * Rz(alpha).Ry(beta).Rz(gamma);
	 */

	R[3 * x + x] =  cos(alpha) * cos(beta) * cos(gamma) - sin(alpha) * sin(gamma);
	R[3 * x + y] = -cos(alpha) * cos(beta) * sin(gamma) - sin(alpha) * cos(gamma);
	R[3 * x + z] =  cos(alpha) * sin(beta);
	R[3 * y + x] =  sin(alpha) * cos(beta) * cos(gamma) + cos(alpha) * sin(gamma);
	R[3 * y + y] = -sin(alpha) * cos(beta) * sin(gamma) + cos(alpha) * cos(gamma);
	R[3 * y + z] =  sin(alpha) * sin(beta);
	R[3 * z + x] = -sin(beta) * cos(gamma);
	R[3 * z + y] =  sin(beta) * sin(gamma);
	R[3 * z + z] =  cos(beta);

	return R;
}


double *euler_inv_rotation_matrix(double gamma, double beta, double alpha)
{
	return euler_rotation_matrix(-alpha, -beta, -gamma);
}


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

	return euler_rotation_matrix(-SKY_MW_Z_J2000_RA_RAD, SKY_MW_Z_J2000_DEC_RAD - M_PI_2, asin(cos(SKY_MW_X_J2000_DEC_RAD) * sin(SKY_MW_Z_J2000_RA_RAD - SKY_MW_X_J2000_RA_RAD)));
}


double *sky_galactic_to_equatorial_rot_matrix(void)
{
	return euler_inv_rotation_matrix(-SKY_MW_Z_J2000_RA_RAD, SKY_MW_Z_J2000_DEC_RAD - M_PI_2, asin(cos(SKY_MW_X_J2000_DEC_RAD) * sin(SKY_MW_Z_J2000_RA_RAD - SKY_MW_X_J2000_RA_RAD)));
}
