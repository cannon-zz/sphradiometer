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
%include carrays.i
%include ccomplex.i
%include cdata.i
%include cpointer.i
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

%inline %{
        struct correlator_plan_fd *pick_ith_correlator_plan_fd(struct correlator_plan_fd **array, int i) {
                return array[i];
        }


        double pick_deltaT_from_COMPLEX16TimeSeries(COMPLEX16TimeSeries *series) {
                return series->deltaT;
        }


        long pick_length_from_COMPLEX16TimeSeries(COMPLEX16TimeSeries *series) {
                return series->data->length;
        }


        void free_SNRTimeSeries(COMPLEX16TimeSeries *series) {
                XLALDestroyCOMPLEX16TimeSeries(series);
        }


        void free_SNRSequence(COMPLEX16Sequence *series) {
                XLALDestroyCOMPLEX16Sequence(series);
        }


        int create_sph_COMPLEX16TimeSeries(COMPLEX16TimeSeries *result, char *name, double epoch_Seconds, double epoch_NanoSeconds, double f0, double deltaT, LALUnit *sampleUnits, size_t length, double complex *data) {
                int i;
                if(!result || !data) {
                        fprintf(stderr, "memory error\n");
                        return -1;
                }
                result->data = XLALCreateCOMPLEX16Sequence(length);

                result->epoch.gpsSeconds = epoch_Seconds;
                result->epoch.gpsNanoSeconds = epoch_NanoSeconds;

                strncpy(result->name, name, LALNameLength - 1);
                result->name[LALNameLength - 1] = '\0';
                result->f0 = f0;
                result->deltaT = deltaT;
                for(i = 0; i < (int) result->data->length; i++) {
                        result->data->data[i] = data[i];
                }
                if(sampleUnits)
                        result->sampleUnits = *sampleUnits;

                return 0;
        }


        int create_sph_COMPLEX16Sequence(COMPLEX16Sequence *result, size_t length, double complex *data) {
                int i;
                if(!result || !data) {
                        fprintf(stderr, "memory error\n");
                        return -1;
                }

                result->data = malloc(length * sizeof(*result->data));
                if(!result->data) {
                        fprintf(stderr, "memory error\n");
                        return -1;
                }

                result->length = length;
                for(i = 0; i < (int) length; i++) {
                        result->data[i] = data[i];
                }

                return 0;
        }
%}

%include <sphradiometer/instrument.h>
%include <sphradiometer/sh_series.h>
%include <sphradiometer/inject.h>
%include <sphradiometer/projection.h>
%include <sphradiometer/correlator.h>
%include <sphradiometer/sky.h>
%include <sphradiometer/deconvolution.h>
%include <sphradiometer/inspiral_sky_map.h>

%pointer_functions(unsigned int, uintp);
%pointer_functions(struct correlator_network_plan_fd, correlator_network_plan_fdp);
%pointer_functions(struct autocorrelator_network_plan_fd, autocorrelator_network_plan_fdp);
%pointer_functions(struct sh_series, sh_seriesp);
%pointer_functions(struct sh_series *, sh_seriespp);
%pointer_functions(COMPLEX16TimeSeries, COMPLEX16TimeSeriesp);
%pointer_functions(COMPLEX16Sequence, COMPLEX16Sequencep);

%array_functions(double, double_array);
%array_functions(double *, doublep_array);
%array_functions(double complex, double_complex_array);
%array_functions(COMPLEX16TimeSeries *, COMPLEX16TimeSeries_array);
%array_functions(COMPLEX16Sequence *, COMPLEX16Sequence_array);
