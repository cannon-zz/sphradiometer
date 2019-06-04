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


#include <complex>
#include <unistd.h>
#include <alm.h>
#include <alm_fitsio.h>
#include <fitshandle.h>
#include <sphradiometer/sh_series.h>


/*
 * ============================================================================
 *
 *                              Type Conversion
 *
 * ============================================================================
 */


/* NOTE:  the coefficients are actually layed out in memory identically, so
 * this could be done more efficiently with a direct memory copy if I knew
 * C++ */


/* NOTE:  although I cannot find it stated anywhere in the healpix
 * documentation, I believe the entire library only supports real-valued
 * functions on the sphere.  the library does not seem to record the
 * co-efficients for m < 0, attempting to set them causes memory corruption
 * related errors.  it would be good to add some error checking to confirm
 * that the sh_series objects being written to disk obey the symmetries
 * required by real-valued functions on the sphere, but it's tricky:  what
 * tolerance do you demand?  for now, it's left as an exercise for the user
 * to remember this limitation of healpix. */


static Alm<xcomplex<double>> *sh_series_to_healpix_Alm(const struct sh_series *series)
{
	int m_max = series->polar ? 0 : series->l_max;
	Alm<xcomplex<double>> *alm = new Alm<xcomplex<double>>(series->l_max, m_max);
	int l, m;

	/* FIXME:  m < 0 is messed up */
	for(m = 0; m <= +m_max; m++)
		for(l = abs(m); l <= (int) series->l_max; l++)
			(*alm)(l, m) = sh_series_get(series, l, m);

	return alm;
}


static struct sh_series *sh_series_from_healpix_Alm(Alm<xcomplex<double>> &alm)
{
	struct sh_series *series;
	int l_max = alm.Lmax();
	int m_max = alm.Mmax();
	int l, m;

	if(m_max == 0)
		series = sh_series_new(l_max, 1);
	else if(m_max == l_max)
		series = sh_series_new(l_max, 0);
	else
		return NULL;
	if(!series)
		return NULL;

	/* FIXME:  m < 0 is messed up */
	for(m = 0; m <= +m_max; m++)
		for(l = abs(m); l <= (int) series->l_max; l++) {
			std::complex<double> x = alm(l, m);
			sh_series_set(series, l, m, *(double _Complex *) &x);
			if(m) {
				x = conj(x) * ((m & 1) ? -1. : +1.);
				sh_series_set(series, l, -m, *(double _Complex *) &x);
			}
		}

	return series;
}


/*
 * ============================================================================
 *
 *                              Type Conversion
 *
 * ============================================================================
 */


extern "C"
int sh_series_write_healpix_alm(const struct sh_series *series, const char *filename)
{
	Alm<xcomplex<double>> *alms;
	fitshandle f;

	try {
		alms = sh_series_to_healpix_Alm(series);
	} catch (std::exception &e) {
		return -1;
	}

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

	try {
		f.create(filename);
		write_Alm_to_fits(f, *alms, alms->Lmax(), alms->Mmax(), planckType<double>());
	} catch (std::exception &e) {
		perror(e.what());
		delete alms;
		return -1;
	}
	delete alms;

	return 0;
}


extern "C"
struct sh_series *sh_series_read_healpix_alm(const char *filename)
{
	fitshandle f;
	int l_max;
	int m_max;
	Alm<xcomplex<double>> *alms = NULL;
	struct sh_series *series;

	try {
		f.open(filename);
		get_almsize(f, l_max, m_max);
		if(l_max < 1 || (m_max != 0 && m_max != l_max))
			throw std::runtime_error("unsupported (l_max, m_max)");
		alms = new Alm<xcomplex<double>>(l_max, m_max);
		read_Alm_from_fits(f, *alms, l_max, m_max);
	} catch (std::exception &e) {
		perror(e.what());
		delete alms;
		return NULL;
	}

	series = sh_series_from_healpix_Alm(*alms);

	delete alms;

	return series;
}
