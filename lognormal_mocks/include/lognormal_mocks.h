#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <symmetric_legendre.h>
#include <healpix_cxx_help.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void lognormal_mocks_stats_fullsky(int Nmap,
				   double *rhobar, // mean density (input)
				   int Nl, double *Cldelta, // power spectrum of overdensity (input)
				   double *gaussbar, double *Clgauss,  // mean and overdensity of Gaussianized field (output)
				   double *xidelta, int Nth // correlation function of overdensity
				   );

void lognormal_mocks_cholesky_spectra(int Nmaps,
				      int Nl, double *Cl,
				      double *L);
void lognormal_mocks_invert_spectra(int Nmaps,
				    int Nl, double *L,
				    double *Sinv);

void lognormal_mocks_invert_spectra_lorder(int Nmaps,
				    int Nl, double *S,
					   double *Sinv);

void lognormal_mocks_genalm(gsl_rng *r, int Nmaps, int Nl, double *L, complex *alm);
void lognormal_mocks_synthesize_fullsky(complex *almset, int Nmaps, int Nl, int nside, double *mapset);
void lognormal_mocks_mapgauss2lognormal(double *gaussbar, double *mapgauss, int Nmaps, int nside, double *maplognormal);
void lognormal_mocks_density2poisson(double *map, gsl_rng *r, int Nmaps, int nside, double *mappoisson);

// Lapack cholesky routine

void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO );
void dpotri_(char *UPLO, int *N, double *A, int *LDA, int *INFO );

