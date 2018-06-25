/* This code takes means and cross-power spectra as inputs, then outputs means and cross-spectra for the Gaussianized fields gauss0 = log(map0), gauss1 = log(map1) ... */

#include <lognormal_mocks.h>


void lognormal_mocks_stats_fullsky(int Nmaps,
				   double *rhobar, // mean density (input)
				   int Nl, double *Cl, // power spectrum of map (input)
				   double *gaussbar, double *Clgauss,  // mean and overdensity of Gaussianized field (output)
				   double *xidelta, int Nth // correlation function of overdensity
				   ) {
  /* Nth = Num of theta controls the accuracy of the intermediate correlations functions (xidelta & xigauss). */

  double vargauss;
  double *xigauss = calloc(Nmaps*Nmaps*Nth,sizeof(double)); // Correlation function of Gaussianized field
  double *Cldelta = calloc(Nmaps*Nmaps*Nl,sizeof(double)); // Power spectrum of Gaussianized field

  int il;
  int ith;

  // make theta and ell arrays;
  double *theta=calloc(Nth,sizeof(double));
  double *ell=calloc(Nl,sizeof(double));
  for (il=0;il<Nl; il++) {
    ell[il] = il;
  }
  for (ith=0;ith<Nth; ith++) {
    theta[ith] = M_PI/(Nth-1)*ith;
  }

  
  
  int B,A;
 
  for (B=0;B<Nmaps;B++) {
    for (A=0;A<=B;A++) {
      printf("Gaussianizing %d %d\n",B,A);

      
      // Compute power spectrum of overdensity
      for (il=0;il<Nl; il++) {
	Cldelta[B*Nmaps*Nl + A*Nl + il] = Cl[B*Nmaps*Nl + A*Nl + il]/rhobar[B]/rhobar[A];
      }
      
      // Compute correlation function of overdensity
      Cl2W_nospline(theta, &xidelta[B*Nmaps*Nth + A*Nth], Nth, &Cldelta[B*Nmaps*Nl + A*Nl], Nl);

      // Compute cross correlations of Gaussianized field
      for (ith=0;ith<Nth; ith++) {
	xigauss[B*Nmaps*Nth + A*Nth + ith] = log( 1 + xidelta[B*Nmaps*Nth + A*Nth + ith] );
      }

      /*fp = fopen("output/test_lognormal_mocks/xigauss.dat","w");
      for (ith=0;ith<Nth;ith++) {
	fprintf(fp,"%e %e %e %e\n",M_PI/(Nth-1)*ith,
		xigauss[0*Nmaps*Nth + 0*Nth + ith],
		xigauss[1*Nmaps*Nth + 0*Nth + ith],
		xigauss[1*Nmaps*Nth + 1*Nth + ith]);
      }
      fclose(fp);*/


      // Compute power spectrum of Gaussianized overdensity
      W2Cl(theta, &xigauss[B*Nmaps*Nth + A*Nth], Nth, ell, &Clgauss[B*Nmaps*Nl + A*Nl], Nl);	

      // Compute mean of Gaussianized field
      if (A==B) {
	vargauss = xigauss[B*Nmaps*Nth + A*Nth + 0];
	gaussbar[B] = log(rhobar[B]) - 0.5*vargauss; 
      }
    }
  }


  // Rescaling Gaussianized power spectrum not necessary.


  free(xigauss);
  free(Cldelta);
  free(theta);

}


void lognormal_mocks_cholesky_spectra(int Nmaps,
				      int Nl, double *Cl,
				      double *L) {
  
  // Generate the Cholesky lower triangles for the correlated spectra Cl
  // Only uses the lower triangle of Cl
  //
  // Note that the ordering of Cl and L are different: Cl keep alike spectra together, L keeps alike ells together.
  
  double *Ltemp = calloc(Nmaps*Nmaps, sizeof(double));
  
  int A,B,il;
  int info;
  char uplo = 'U';
  
  for (il=0;il<Nl; il++) {

    for (B=0;B<Nmaps;B++) {
      for (A=0;A<=B;A++) {
	Ltemp[B*Nmaps + A] = Ltemp[A*Nmaps + B] = Cl[B*Nmaps*Nl + A*Nl + il];
      }
    }

    
    // Cholesky decompose with lapack
    dpotrf_( &uplo, &Nmaps, Ltemp, &Nmaps, &info );
    
    if (info == 0) {
      for (B=0;B<Nmaps;B++) {
	for (A=0;A<=B;A++) {
	  L[il*Nmaps*Nmaps + B*Nmaps + A] = Ltemp[B*Nmaps + A];
	}
      }
    } else {
      for (B=0;B<Nmaps;B++) {
	for (A=0;A<=B;A++) {
	  L[il*Nmaps*Nmaps + B*Nmaps + A] = 0.0;
	}
      }
    }
    
    
  }

  free(Ltemp);
}


void lognormal_mocks_invert_spectra(int Nmaps,
				    int Nl, double *L,
				    double *Sinv) {
  // Invert from the Cholesky decomposition
  // Like the Cholesky, the inverse keeps like ells together
  
  double *Sinvtemp = calloc(Nmaps*Nmaps, sizeof(double));
  int A,B,il;
  int info;
  char uplo = 'U';

  for (il=0;il<Nl; il++) {

    for (B=0;B<Nmaps;B++) {
      for (A=0;A<=B;A++) {
	Sinvtemp[B*Nmaps + A] = Sinvtemp[A*Nmaps + B] = L[il*Nmaps*Nmaps + B*Nmaps + A];
      }
    }

    dpotri_( &uplo, &Nmaps, Sinvtemp, &Nmaps, &info );

    
    if (info == 0) {
      for (B=0;B<Nmaps;B++) {
	for (A=0;A<=B;A++) {
	  Sinv[il*Nmaps*Nmaps + A*Nmaps + B] = Sinv[il*Nmaps*Nmaps + B*Nmaps + A] = Sinvtemp[B*Nmaps + A];
	}
      }
    } else {
      for (B=0;B<Nmaps;B++) {
	for (A=0;A<=B;A++) {
	  Sinv[il*Nmaps*Nmaps + A*Nmaps + B] = Sinv[il*Nmaps*Nmaps + B*Nmaps + A] = 0.0;
	}
      }
    }
    
  }
    

  free(Sinvtemp);

}

void lognormal_mocks_invert_spectra_lorder(int Nmaps,
				    int Nl, double *S,
				    double *Sinv) {
  // Invert a whole set of cross spectrum
  // in this version, both initial spectrum and the inverse keeps like ells together
  
  double *Sinvtemp = calloc(Nmaps*Nmaps, sizeof(double));
  int A,B,il;
  int info;
  char uplo = 'U';

  for (il=0;il<Nl; il++) {

    for (B=0;B<Nmaps;B++) {
      for (A=0;A<=B;A++) {
	Sinvtemp[B*Nmaps + A] = Sinvtemp[A*Nmaps + B] = S[il*Nmaps*Nmaps + B*Nmaps + A];
      }
    }

    dpotrf_( &uplo, &Nmaps, Sinvtemp, &Nmaps, &info );
    dpotri_( &uplo, &Nmaps, Sinvtemp, &Nmaps, &info );

    
    if (info == 0) {
      for (B=0;B<Nmaps;B++) {
	for (A=0;A<=B;A++) {
	  Sinv[il*Nmaps*Nmaps + A*Nmaps + B] = Sinv[il*Nmaps*Nmaps + B*Nmaps + A] = Sinvtemp[B*Nmaps + A];
	}
      }
    } else {
      for (B=0;B<Nmaps;B++) {
	for (A=0;A<=B;A++) {
	  Sinv[il*Nmaps*Nmaps + A*Nmaps + B] = Sinv[il*Nmaps*Nmaps + B*Nmaps + A] = 0.0;
	}
      }
    }
    
  }
    

  free(Sinvtemp);

}



void lognormal_mocks_genalm(gsl_rng *r, int Nmaps, int Nl, double *L, complex *alm){

  int l,m,i;
  complex *u = calloc(Nmaps,sizeof(complex)); 
  int A,B;
  int Nalm = Nl*Nl;
  
  for (l=0;l<Nl;l++) {
 
    for (m=0; m<=l; m++) {
      i = l*l + l + m;

      // Make uniform deviates
      if (m==0) {
	for (A=0;A<Nmaps;A++) {
	  u[A] = gsl_ran_gaussian_ziggurat(r,1.0);
	}
      } else {
	for (A=0;A<Nmaps;A++) {
	  u[A] = gsl_ran_gaussian_ziggurat(r,M_SQRT1_2) + I*gsl_ran_gaussian_ziggurat(r,M_SQRT1_2);
	}
      }

      // project them into the alms with matrix multiplication
      for (B=0;B<Nmaps;B++) {
	alm[B*Nalm + i] = 0.0;
	for (A=0;A<=B;A++) { // lower triangle
	  alm[B*Nalm + i] += L[l*Nmaps*Nmaps + B*Nmaps + A] * u[A];
	}
      }
      
    }
  }


  free(u);
}

void lognormal_mocks_synthesize_fullsky(complex *almset, int Nmaps, int Nl, int nside, double *mapset) {
  int B;
  char filename[1024];
  int lmax = Nl - 1;
  int Nalm = Nl*Nl;
  int npix = 12*nside*nside;
  
  for (B=0;B<Nmaps;B++) {
    hpcxx_alm2map(&almset[B*Nalm], lmax, lmax, nside, &mapset[B*npix]);
  }


}


void lognormal_mocks_mapgauss2lognormal(double *gaussbar, double *mapgauss, int Nmaps, int nside, double *maplognormal) {

  int B;
  int p,npix = 12*nside*nside;

  
  for (B=0;B<Nmaps;B++) {
    for (p=0; p<npix;p++) {
      maplognormal[B*npix + p] = exp ( gaussbar[B] + mapgauss[B*npix + p] );
    }
  }
  
}


void lognormal_mocks_density2poisson(double *map, gsl_rng *r, int Nmaps, int nside, double *mappoisson){

  int B;
  int p,npix = 12*nside*nside;
  double Omegapix = 4*M_PI/npix;
  double mu;
  
  for (B=0;B<Nmaps;B++) {
    for (p=0; p<npix;p++) {
      mu = map[B*npix + p] * Omegapix;
      mappoisson[B*npix + p] = gsl_ran_poisson (r, mu);
    }
  }
  



}
