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


#ifndef __RADIOMETER_SKY_H__
#define __RADIOMETER_SKY_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


/*
 * ============================================================================
 *
 *                                   Macros
 *
 * ============================================================================
 */


/*
 * Equatorial co-ordinates for the north pole ("z axis") and co-ordinate
 * origin ("x axis") of the galactic co-ordinate system.
 *
 * Reference:
 *
 * 	Wikipedia (https://en.wikipedia.org/wiki/Galactic_coordinate_system)
 */


#define SKY_MW_Z_J2000_RA_RAD 3.3658674624710647 /* = 12h 51.4m */
#define SKY_MW_Z_J2000_DEC_RAD 0.4735078260660616 /* = 27.13o */
#define SKY_MW_X_J2000_RA_RAD 4.649557127312894 /* = 17 h 45.6 m */
#define SKY_MW_X_J2000_DEC_RAD -0.5050982855271591 /* = -28.94o */


/*
 * Equatorial co-ordinates for the galactic centre.  NOTE:  Radio source
 * Sagittarius A* is not at these co-ordinates, but it is within 0.1o of
 * these co-ordinates which is the error budget for the 1958 co-ordinate
 * system's definition
 *
 * Reference:
 *
 * 	Wikipedia (https://en.wikipedia.org/wiki/Galactic_coordinate_system)
 */


#define SKY_MW_CENTRE_J2000_RA_RAD 4.649557127312894 /* = 17 h 45.6 m */
#define SKY_MW_CENTRE_J2000_DEC_RAD -0.5050982855271591 /* = -28.94o */


/*
 * Equatorial co-ordinates for the Crab pulsar.
 *
 * Reference:
 *
 * 	ATNF pulsar database (http://www.atnf.csiro.au/research/pulsar/psrcat)
 */


#define SKY_CRAB_PULSAR_J2000_RA_RAD 1.4596750675891823 /* = 5h 34m 31.973 s */
#define SKY_CRAB_PULSAR_J2000_DEC_RAD 0.38422482944113812 /* = +22o 00' 52.06" */


/*
 * ============================================================================
 *
 *                                 Prototypes
 *
 * ============================================================================
 */


double *euler_rotation_matrix(double, double, double);
double *euler_inv_rotation_matrix(double, double, double);


double *sky_equatorial_to_galactic_rot_matrix(void);
double *sky_galactic_to_equatorial_rot_matrix(void);


#endif  /* __RADIOMETER_SKY_H__ */
