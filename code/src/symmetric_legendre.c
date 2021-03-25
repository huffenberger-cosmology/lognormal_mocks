/*
    This code is part of lognormal_mocks.
    Copyright (C) 2021 Kevin M. Huffenberger, khuffenberger@fsu.edu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/



#include <symmetric_legendre.h>


void W2Cl(double *theta, double *W, int Ntheta, double *l, double *C, int Nl) {
  // Take legendre transform w/ Riemann sum

  int lmax = l[Nl-1];
  int mmax = 0;
  int il,ith;
 
  hpl_init_rescale(); 
  double *mfac = hpl_alloc_mfac(mmax);
  hpl_gen_mfac(mmax, mfac);
  double *recfac = hpl_alloc_recfac(lmax);
  hpl_gen_recfac(lmax, mmax, recfac);
  
  double *lam_l0 = calloc((lmax+1)*Ntheta,sizeof(double));
  
  //  Generate lam_l0 = sqrt((2l+1)/4pi) P_l (cos(theta))
  for (ith = 0; ith<Ntheta;ith++) {
    double cth = cos(theta[ith]);
    double sth = sin(theta[ith]);
    
    hpl_do_lam_lm(lmax, mmax, cth, sth, mfac[mmax], recfac, &lam_l0[ith*(lmax+1)]);
  }

  double norm = 2*M_PI;
  for (il=0;il<Nl;il++) {
    double Pfac = sqrt((2.0*l[il]+1.0)/4.0/M_PI);
    C[il] = 0;

    //   printf("");

    for (ith = 0; ith<Ntheta-1;ith++) {
      double dtheta =  theta[ith+1] - theta[ith];
      double Pl = lam_l0[ith*(lmax+1) + (int)l[il]] / Pfac;
      
      if (theta[ith]==0.0) Pl = 1.0;

      //   if (il==3) printf("l %e dtheta %e Pl %e\n",l[il],dtheta,Pl);

      C[il] += dtheta * sin(theta[ith]) * W[ith] * Pl;
    }

    C[il] *= norm;
  }


  free(mfac);
  free(recfac);
  free(lam_l0);
}



void Cl2W(double *theta, double *W, int Ntheta, double *l, double *C, int Nl) {
  // Take legendre transform.  Sparsely sample Cl's are splined to give value at every l.

  gsl_spline *Cspline = gsl_spline_alloc (gsl_interp_cspline, Nl);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline_init(Cspline, l, C, Nl);

  int lmax = l[Nl-1];
  int il;

  double *Call = calloc(lmax+1,sizeof(double));
  for (il=0;il<=lmax;il++) {
    Call[il] = gsl_spline_eval(Cspline, il, acc);
  }
  
  Cl2W_nospline(theta, W, Ntheta, Call, lmax+1);

  free(Call);
  gsl_spline_free (Cspline);
  gsl_interp_accel_free(acc);
}

void Cl2W_nospline(double *theta, double *W, int Ntheta, double *C, int Nl) {

  int lmax = Nl - 1;
  int mmax = 0;
  int ith;
  int il;

  
  hpl_init_rescale(); 
  double *mfac = hpl_alloc_mfac(mmax);
  hpl_gen_mfac(mmax, mfac);
  double *recfac = hpl_alloc_recfac(lmax);
  hpl_gen_recfac(lmax, mmax, recfac);
  
  double *lam_l0 = calloc((lmax+1)*Ntheta,sizeof(double));
  
  //  Generate lam_l0 = sqrt((2l+1)/4pi) P_l (cos(theta))
  for (ith = 0; ith<Ntheta;ith++) {
    double cth = cos(theta[ith]);
    double sth = sin(theta[ith]);
    
    hpl_do_lam_lm(lmax, mmax, cth, sth, mfac[mmax], recfac, &lam_l0[ith*(lmax+1)]);
  }





  double norm = 1.0;
  for (ith = 0; ith<Ntheta-1;ith++) {
    W[ith] = 0;

    for (il=0;il<=lmax;il++) {
      double Pfac = sqrt((2.0*il+1.0)/4.0/M_PI);
      
      //   printf("");

      double Pl = lam_l0[ith*(lmax+1) + il] / Pfac;
      if (theta[ith]==0.0) Pl = 1.0;

      double Cl = C[il];

      W[ith] += Cl * Pfac*Pfac * Pl;
	
    }
    W[ith] *= norm;

  }

  free(mfac);
  free(recfac);
  free(lam_l0);

}


  
// To interpolate W or Cl 
/* gsl_spline *Wspline = gsl_spline_alloc (gsl_interp_cspline, Ntheta); */
/* gsl_interp_accel *acc= gsl_interp_accel_alloc(); */
/* gsl_spline_init(Wspline, theta, W, Ntheta); */
/* //  gsl_spline_eval(Wspline, th, acc); */
/* gsl_spline_free (Wspline); */
