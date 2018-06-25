#include <healpix_legendre_c.h>

double *hpl_alloc_recfac(int l_max){
  // From gen_recfac:
  //   real(DP),     intent(OUT), dimension(0:1, 0:l_max)  :: recfac

  return( (double *) calloc(2*(l_max+1),sizeof(double)) );
}

double *hpl_alloc_mfac(int m_max){
  // From gen_mfac:
  //   real(DP),     intent(OUT), dimension(0:m_max)  :: m_fact

  return( (double *) calloc(m_max+1,sizeof(double)) );
}


void hpl_do_lam_lm(int lmax, int m, double cth, double sth, double mfac, double *recfac, double *lam_lm) {
  healpix_do_lam_lm(&lmax, &m, &cth, &sth, &mfac, recfac, lam_lm);
}

void hpl_gen_mfac(int m_max, double *m_fact) {
  healpix_gen_mfac(&m_max, m_fact);
}

void hpl_gen_recfac(int l_max, int m, double *recfac) {
  healpix_gen_recfac(&l_max, &m, recfac);
}
