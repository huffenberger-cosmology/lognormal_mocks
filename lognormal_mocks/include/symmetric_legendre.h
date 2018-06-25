#include <healpix_legendre_c.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>


void W2Cl(double *theta, double *W, int Ntheta, double *l, double *C, int Nl);
void Cl2W(double *theta, double *W, int Ntheta, double *l, double *C, int Nl);
void Cl2W_nospline(double *theta, double *W, int Ntheta, double *C, int Nl); // uses every ell
