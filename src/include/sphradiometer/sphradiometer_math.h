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
 * This file contains arithmetic macros for use in the library code, and is
 * not meant to provide an exported interface.  This file is *not*
 * installed.
 */


#ifndef __RADIOMETER_MATH_H__
#define __RADIOMETER_MATH_H__


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <complex.h>
#include <math.h>


#ifdef __cplusplus
extern "C" {
#define complex _Complex
#endif


/*
 * ============================================================================
 *
 *                                 Utilities
 *
 * ============================================================================
 */


/*
 * fma() for double complex * double complex + double complex.
 */


#ifndef SWIG
static double complex ccma(double complex, double complex, double complex) __attribute__ ((unused));
#endif
static double complex ccma(double complex x, double complex y, double complex z)
{
	return fma(-cimag(x), cimag(y), fma(creal(x), creal(y), creal(z))) + I * fma(cimag(x), creal(y), fma(creal(x), cimag(y), cimag(z)));
}


/*
 * fma() for double complex * double + double complex.
 */


#ifndef SWIG
static double complex cma(double complex, double, double complex) __attribute__ ((unused));
#endif
static double complex cma(double complex x, double y, double complex z)
{
	return fma(creal(x), y, creal(z)) + I * fma(cimag(x), y, cimag(z));
}


#ifdef __cplusplus
#undef complex
}
#endif


#endif	/* __RADIOMETER_MATH_H__ */
