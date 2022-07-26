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
#include <errno.h>
#include <math.h>
#include <stdexcept>
#include <stdlib.h>
#include <unistd.h>
#include <alm.h>
#include <alm_fitsio.h>
#include <alm_healpix_tools.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
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
 * tolerance do you demand?  for now, we use sh_series_real() to obtain the
 * coefficients of the real-valued part of the function, and it's left as
 * an exercise for the user to remember this limitation of healpix.  */


static Alm< xcomplex<double> > *sh_series_to_healpix_Alm(const struct sh_series *series)
{
	sh_series *real = sh_series_real(sh_series_copy(series));
	int m_max = series->polar ? 0 : series->l_max;
	Alm< xcomplex<double> > *alm = new Alm< xcomplex<double> >(series->l_max, m_max);
	int l, m;

	for(m = 0; m <= +m_max; m++)
		for(l = abs(m); l <= (int) series->l_max; l++)
			(*alm)(l, m) = sh_series_get(real, l, m);

	sh_series_free(real);

	return alm;
}


static struct sh_series *sh_series_from_healpix_Alm(Alm< xcomplex<double> > &alm)
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

	/* NOTE:  loop is only over non-negative m because m < 0 is not
	 * available in Alm objects. */
	for(m = 0; m <= +m_max; m++)
		for(l = abs(m); l <= (int) series->l_max; l++) {
			std::complex<double> x = alm(l, m);
			sh_series_set(series, l, m, *(double _Complex *) &x);
			/* fill in the missing coefficients by assuming a
			 * real-valued function on sphere. */
			if(m) {
				x = (m & 1) ? -conj(x) : conj(x);
				sh_series_set(series, l, -m, *(double _Complex *) &x);
			}
		}

	return series;
}


/*
 * ============================================================================
 *
 *                                    I/O
 *
 * ============================================================================
 */


extern "C"
int sh_series_write_healpix_alm(const struct sh_series *series, const char *filename)
{
	Alm< xcomplex<double> > *alms;
	fitshandle f;

	try {
		alms = sh_series_to_healpix_Alm(series);
	} catch (std::exception &e) {
		return -1;
	}

	/* FITS barfs if the file exists, which is annoying, so delete it
	 * first because people expect functions like this to overwrite the
	 * target file if it exists.  ignore errors from unlink() because
	 * it will complain if the file doesn't exist (it shouldn't exist,
	 * that is expected), and if the delete fails because of some
	 * filesystem malfunction then write_healpfix() will also fail and
	 * we'll let it produce the error message */

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
	Alm< xcomplex<double> > alms;

	try {
		f.open(filename);
		f.goto_hdu(2);
		get_almsize(f, l_max, m_max);
		if(l_max < 0 || (m_max != 0 && m_max != l_max))
			throw std::runtime_error("unsupported (l_max, m_max)");
		read_Alm_from_fits(f, alms, l_max, m_max);
	} catch (std::exception &e) {
		perror(e.what());
		return NULL;
	}

	return sh_series_from_healpix_Alm(alms);
}


/*
 * return the smallest power of 2 not smaller than x.  returns -1 on
 * overflow.
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
 * compute the healpix nside parameter corresponding to an l_max.  the
 * conversion begins by imposing that the healpix map have the correct
 * number of rings of isolatitude for the given l_max, and then rounds up
 * the result to a power 2 because it seems healpix prefers that.  returns
 * -1 on overflow.
 */


static int lmax2nside(int l_max)
{
	/* healpix map has (4 * nside - 1) rings of isolatitude
	 *
	 * from the harmonic transform equations, there must be 2 * (l_max
	 * + 1) rings of isolatitude
	 *
	 * therefore
	 *
	 * 	nside = (2 * l_max + 3) / 4
	 *
	 * but healpix wants nside to be a power of 2, so we round up to
	 * one.  NOTE:  recalling that l_max + 1 is often called the
	 * bandwidth, the relationship between l_max and nside given above
	 * is equivalent to nside = (2 * bandwidth + 1) / 4.
	 */

	return ceilpow2(ceil(l_max / 2. + 0.75));
}


/*
 * Write an sh_series to a real-valued healpix map FITS file.  NOTE:  this
 * is expensive because it includes the cost of a frequency domain -to-
 * pixel domain conversion.
 */


extern "C"
int sh_series_write_healpix_map(const struct sh_series *series, const char *filename)
{
	Healpix_Map< double > map(Healpix_Base::nside2order(lmax2nside(series->l_max)), RING);
	fitshandle f;

	try {
		Alm< xcomplex<double> > *alms = sh_series_to_healpix_Alm(series);
		alm2map(*alms, map);
		delete alms;
	} catch (std::exception &e) {
		return -1;
	}

	/* FITS barfs if the file exists, which is annoying, so delete it
	 * first because people expect functions like this to overwrite the
	 * target file if it exists.  ignore errors from unlink() because
	 * it will complain if the file doesn't exist, which we don't care
	 * about (it shouldn't exist, that is expected), and if the file
	 * does exist and the delete fails because of some filesystem
	 * malfunction then write_healpfix() will also fail and we'll let
	 * it produce some kind of error message */

	{
	int errsv = errno;
	unlink(filename);
	errno = errsv;
	}

	try {
		f.create(filename);
		write_Healpix_map_to_fits(f, map, planckType<double>());
	} catch (std::exception &e) {
		perror(e.what());
		return -1;
	}

	return 0;
}
