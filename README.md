sphradiometer --- Library for synthesis imaging on the sphere
=============================================================

This contains the code used for Cannon, K., "Efficient algorithm for computing the time-resolved full-sky cross power in an interferometer with omnidirectional elements", Phys. Rev. D 75, 123003, 2007-06-13 (https://doi.org/10.1103/PhysRevD.75.123003), and also Tsutsui, T., et al., "High speed source localization in searches for gravitational waves from compact object collisions", Phys. Rev. D 103, 043011, 2021-02-22 (https://doi.org/10.1103/PhysRevD.103.043011).

The first paper generalizes the van Cittert-Zernike theorem away from the narrow-band, small field of view, limit, presenting a numerical approximation scheme for a specific integral that appears in the problem of omnidirectional aperture synthesis in a low-frequency broad-band interferometer.  The problem being considered was gravitational-wave sky mapping, which does not require the sky to be mapped beyond l=25 or so.  In this regime, the algorithm is extraordinarily fast, producing sky maps in tens of milliseconds on a single core CPU.

The second paper applies the algorithm to the specific case of localizing the source of gravitational-waves from a compact object collision, by approximating the solution to the problem making use of the integral form for which the first paper provided a fast numerical scheme.


Installation
============

It is recommend that users install pre-built packages for their system.  .deb and .rpm packages might be available if you look around for them, or you can make them yourself from source.

From git
--------

A git clone requires

	$ ./00init.sh

to be run to generate the build and install framework.  Doing this requires the GNU autotools software suite to be installed in addition to the other build requirements.  After this, follow the instructions for building from source.

From Source
-----------

The software is installed from a source tarball or from a 00init.sh'ed git clone using the standard GNU

	$ ./configure && make install

Expect to provide at least a suitable --prefix option for configure.

Building a Package
------------------

After running configure, a source tarball can be obtained with

	$ make dist

A .deb package set can be obtained with

	$ dpkg-buildpackage --no-pre-clean --no-sign -r fakeroot

If some of the build dependencies have also been installed from source, add --no-check-builddeps to the options, and uncomment the override_dh_shlibdeps target in debian/rules.
