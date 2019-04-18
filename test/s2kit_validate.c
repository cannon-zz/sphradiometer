/*
 * Copyright (C) 2019  Kipp C. Cannon
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


#include <stdio.h>


#include <s2kit/FST_semi_memo.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	int l_max = 10;
	int l, m;

	for(l = 0; l <= l_max; l++)
		for(m = -l; m <= l; m++)
			if(seanindex(m, l, l_max + 1) != (int) sh_series_params_lmoffset(l_max, l, m)) {
				fprintf(stderr, "seanindex(%d, %d, %d) = %d but sh_series_params_lmoffset(%d, %d, %d) = %zu\n", m, l, l_max + 1, seanindex(m, l, l_max + 1), l_max, l, m, sh_series_params_lmoffset(l_max, l, m));
				exit(1);
			}

	return 0;
}
