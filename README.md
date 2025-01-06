# sphradiometer --- Library for synthesis imaging on the sphere

This contains the code used for Cannon, K., "Efficient algorithm for computing the time-resolved full-sky cross power in an interferometer with omnidirectional elements", Phys. Rev. D 75, 123003, 2007-06-13 (https://doi.org/10.1103/PhysRevD.75.123003), and also Tsutsui, T., et al., "High speed source localization in searches for gravitational waves from compact object collisions", Phys. Rev. D 103, 043011, 2021-02-22 (https://doi.org/10.1103/PhysRevD.103.043011).

The first paper generalizes the van Cittert-Zernike theorem away from the narrow-band, small field of view, limit, presenting a numerical approximation scheme for a specific integral that appears in the problem of omnidirectional aperture synthesis in a low-frequency broad-band interferometer.  The problem being considered was gravitational-wave sky mapping, which does not require the sky to be mapped beyond l=25 or so.  In this regime, the algorithm is extraordinarily fast, producing sky maps in tens of milliseconds on a single core CPU.

The second paper applies the algorithm to the specific case of localizing the source of gravitational-waves from a compact object collision, by approximating the solution to the problem making use of the integral form for which the first paper provided a fast numerical scheme.
