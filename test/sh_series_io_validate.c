/*
 * Copyright (C) 2021  Kipp C. Cannon
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

#include <sphradiometer/diagnostics.h>
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
	struct sh_series *in;
	struct sh_series *out;

	out = sh_series_real(random_sh_series(10, 0));

	fprintf(stderr, "write fits file\n");
	if(sh_series_write_healpix_alm(out, "test.fits")) {
		fprintf(stderr, "write \"test.fits\" failed\n");
		sh_series_free(out);
		exit(1);
	}

	fprintf(stderr, "read fits file\n");
	in = sh_series_read_healpix_alm("test.fits");

	if(sh_series_cmp(out, in)) {
		fprintf(stderr, "data read not equal to data written\n");
		sh_series_free(out);
		sh_series_free(in);
		exit(1);
	}

	sh_series_free(out);
	sh_series_free(in);

	exit(0);
}
