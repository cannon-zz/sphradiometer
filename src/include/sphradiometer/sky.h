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
 *                                 Data Types
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
 * 	Wikipedia (http://en.wikipedia.org/wiki/Galactic_latitude)
 */


#define SKY_MW_Z_RA_RAD 3.3660334141941082 /* = 12h 51m 26.282s */
#define SKY_MW_Z_DEC_RAD 0.47347878572656316 /* = 27o 07' 42.01" */
#define SKY_MW_X_J2000_RA_RAD 4.6496461391047452 /* = 17h 45m 37.224s */
#define SKY_MW_X_J2000_DEC_RAD -0.50503152668327012 /* = -28o 56' 10.23" */


/*
 * Equatorial co-ordinates for the galactic centre.
 *
 * Reference:
 *
 * 	Wikipedia (http://en.wikipedia.org/wiki/Galactic_latitude)
 */


#define SKY_MW_CENTRE_J2000_RA_RAD 4.6498509244036468 /* = 17 h 45 m 40.04 s */
#define SKY_MW_CENTRE_J2000_DEC_RAD -0.50628171572274738 /* = -29o 00' 28.1" */


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


double *sky_rotation_matrix(double, double, double);
double *sky_inv_rotation_matrix(double, double, double);


double *sky_equatorial_to_galactic_rot_matrix(void);
double *sky_galactic_to_equatorial_rot_matrix(void);


#endif  /* __RADIOMETER_SKY_H__ */
