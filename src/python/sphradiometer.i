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
        double **flat_psd_array(int n_ifo, int size) {
                int i, j;
                double **res = malloc(size * sizeof(*res));

                if(!res) {
                        fprintf(stderr, "out of memory\n");
                        return NULL;
                }

                for(i = 0; i < size; i++) {
                        res[i] = malloc(n_ifo * sizeof(*res[i]));
                        if(!res[i]) {
                                fprintf(stderr, "out of memory\n");
                                free(res);
                                return NULL;
                        }
                        for(j = 0; j < n_ifo; j++)
                                res[i][j] = 1;
                }

                return res;
        }


        void free_psd_array(double **array, int size) {
                int i;
                for(i = 0; i < size; i++)
                        free(array[i]);
                free(array);
        }


        struct correlator_plan_fd *pick_ith_correlator_plan_fd(struct correlator_plan_fd **array, int i) {
                return array[i];
        }


        double pick_deltaT_from_COMPLEX16TimeSeries(COMPLEX16TimeSeries *series) {
                return series->deltaT;
        }


        long pick_length_from_COMPLEX16TimeSeries(COMPLEX16TimeSeries *series) {
                return series->data->length;
        }


        void free_SNRTimeSeries(COMPLEX16TimeSeries **series, COMPLEX16Sequence **nseries, int ndet) {
                int k;
                for(k = 0; k < ndet; k++) {
                        XLALDestroyCOMPLEX16TimeSeries(series[k]);
                        XLALDestroyCOMPLEX16Sequence(nseries[k]);
                }

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

%array_functions(double, double_array);
%array_functions(double *, doublep_array);
%array_functions(double complex, double_complex_array);
%array_functions(COMPLEX16TimeSeries *, COMPLEX16TimeSeries_array);
%array_functions(COMPLEX16Sequence *, COMPLEX16Sequence_array);
