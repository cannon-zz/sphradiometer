/*
 * Copyright (C) 2006  Kipp C. Cannon
 * Copyright (C) 2020  Takuya Tsutsui
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

%module sphradiometer
%include ccomplex.i
%include typemaps.i

/*%include <lal/SWIGCommon.i>
#ifndef SWIGIMPORTED
%import <lal/swiglal.i>
#endif*/

%{
#define SWIG_FILE_WITH_INIT
#include <sphradiometer/instrument.h>
#include <sphradiometer/sh_series.h>
#include <sphradiometer/inject.h>
#include <sphradiometer/projection.h>
#include <sphradiometer/correlator.h>
#include <sphradiometer/sky.h>
#include <sphradiometer/deconvolution.h>
#include <sphradiometer/inspiral_sky_map.h>
%}

%include <sphradiometer/instrument.h>
%include <sphradiometer/sh_series.h>
%include <sphradiometer/inject.h>
%include <sphradiometer/projection.h>
%include <sphradiometer/correlator.h>
%include <sphradiometer/sky.h>
%include <sphradiometer/deconvolution.h>
%include <sphradiometer/inspiral_sky_map.h>
