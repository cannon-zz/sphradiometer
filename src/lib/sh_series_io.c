/*
 * Copyright (C) 2006--2009,2019  Kipp C. Cannon
 * Copyright (C) 2019  Takuya Tsutsui
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


#include <chealpix.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                                  File I/O
 *
 * ============================================================================
 */


/*
 * Print the non-zero coefficients in an sh_series object.
 */


void sh_series_print(FILE *f, const struct sh_series *series)
{
	int l, m;
	int lines = 0;

	for(l = 0; l <= (int) series->l_max; l++)
		for(m = series->polar ? 0 : -l; m <= (series->polar ? 0 : l); m++) {
			complex double c = sh_series_get(series, l, m);
			if(c != 0.0) {
				fprintf(f, "(%d,%d) = %.17g + I %.17g\n", l, m, creal(c), cimag(c));
				lines++;
			}
		}
	/* if all coefficients are zero, print at least something to
	 * simplify code that tries to parse this output */
	if(!lines)
		fprintf(f, "(%d,%d) = %.17g + I %.17g\n", 0, 0, 0.0, 0.0);
}


/*
 * return the smallest power of 2 not smaller than x
 */


static int ceilpow2(int x)
{
	int i;
	for(i = 1; i > 0; i <<= 1)
		if(i >= x)
			return i;
	return -1;	/* overflow */
}


/*
 * Write an sh_series to a real-valued healpix map FITS file.
 */


int sh_series_write_healpix_map(const struct sh_series *series, const char *filename)
{
	struct sh_series_eval_interp *interp = sh_series_eval_interp_new(series);
	int nside = ceilpow2(ceil(series->l_max / 2.));
	int npix = nside2npix(nside);
	float *map = malloc(npix * sizeof(*map));
	int ipring;

	if(!interp || !map) {
		sh_series_eval_interp_free(interp);
		free(map);
		return -1;
	}

	for(ipring = 0; ipring < npix; ipring++) {
		double theta, phi;
		pix2ang_ring(nside, ipring, &theta, &phi);
		map[ipring] = creal(sh_series_eval_interp(interp, theta, phi));
	}

	sh_series_eval_interp_free(interp);

	/* FITS barfs if the file exists, which is annoying, so delete it
	 * first because people expect functions like this to overwrite the
	 * target file if it exists.  ignore errors from unlink() because
	 * the file might not exist, and if the delete fails then
	 * write_healpfix() will fail and we'll let it produce the error
	 * message */

	{
	int errsv = errno;
	unlink(filename);
	errno = errsv;
	}

	write_healpix_map(map, nside, filename, 0, "C");

	free(map);

	return 0;
}
