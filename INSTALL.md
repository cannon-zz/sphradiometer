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
