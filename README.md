# lognormal_mocks

This code allows for the production of 2-d correlated lognormal random fields.  These have application as e.g. cosmological (rho > 0) density fields projected on the sky.  The basic math behind this method comes from Carron & Neyrinck (2012) ApJ 750 28, section 3.1.

The test/ directory provides several examples of the code in use.  These codes follow a similar structure.  The mean and the power spectrum is specified for the lognormal field and fed to the function lognormal_mocks_stats().  This produces the mean and power spectrum of the "Gaussianized" field.  We use healpy to generate the Gaussianized maps, then exponentiate them to produce the lognormal maps.  Finally we verify that the lognormal field have the desired power spectra.

'basic.py' makes a lognormal random field in a healpix map with a toy power spectrum

'correlated.py' makes three lognormal maps with specified auto and cross spectra.

'mock_gg.py' makes a realization of a galaxy count and overdensity map, and includes Poisson sampling of the galaxies on top of the lognormal galaxy density field.

'mock_kappag.py' makes correlated galaxy overdensity and CMB lensing kappa maps, including the Poisson sampling for galaxies.

Gabriela Marques provided theory C_l curves for the galaxy overdensity and CMB lensing kappa examples.

The accuracy of the method is set by the maximum multipole l of the input power spectrum and the Ntheta parameter.  Internally, lognormal_mocks_stats() does a Legendre transformation to convert the power spectrum to a correlation function.  These parameters describe the bandwidth and sampling included, and determine the accuracy and the runtime.  Because the transformation in the exponential function is nonlinear, the usual intuition about l_max, Ntheta, and Nside may not hold, and you may need to go to higher resolution than you expect.  Systems where the mean is small compared to the power spectrum are the most non-Gaussian and the most challenging.

The code doesn't do much internal checking, so it is vital that you verify the results match what you expect.  If you need help, ask!

## Code and Dependences
The code is written with the Python C API.

The code depends on GSL libraries and Lapack for interpolation and matrix operations (Cholesky decomposition).  It depends on the Healpix Fortran library (libhealpix) to  conduct the Legendre transform.  Healpix in turn depends on cfitsio and in some installations an OpenMP library like libgomp.
