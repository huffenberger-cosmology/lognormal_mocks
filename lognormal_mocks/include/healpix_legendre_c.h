/////////////////////////////////////////
//
// C wrappers for Healpix Library Assoc. Legendre Functions
//
/////////////////////////////////////////
//
//  Sample calling sequence:
//       hpl_init_rescale();  // all values zero without this.
//
//       double *mfac = hpl_alloc_mfac(mmax);
//       hpl_gen_mfac(mmax, mfac);
//
//       double *recfac = hpl_alloc_recfac(lmax);
//       hpl_gen_recfac(lmax, m, recfac);
//
//       hpl_do_lam_lm(lmax, m, cth, sth, mfac[m], recfac, lam_lm);
//
//  hpl_do_lam_lm computes all l's up to lmax for given m, theta.
//  cth = cos(theta), sth = sin(theta)
//  lam_lm[l] contains lambda_lm( cos(theta) )
//


#include <stdio.h>
#include <stdlib.h>


#define healpix_do_lam_lm __alm_tools_MOD_do_lam_lm
#define healpix_gen_mfac __alm_tools_MOD_gen_mfac
#define healpix_gen_recfac __alm_tools_MOD_gen_recfac
#define hpl_init_rescale __alm_tools_MOD_init_rescale


void hpl_init_rescale();

void hpl_do_lam_lm(int lmax, int m, double cth, double sth, double mfac, double *recfac, double *lam_lm);
void healpix_do_lam_lm(int *lmax, int *m, double *cth, double *sth, double *mfac, double *recfac, double *lam_lm);

void hpl_gen_mfac(int m_max, double *m_fact);
void healpix_gen_mfac(int *m_max, double *m_fact);
double *hpl_alloc_mfac(int m_max);

void hpl_gen_recfac(int l_max, int m, double *recfac);
void healpix_gen_recfac(int *l_max, int *m, double *recfac);
double *hpl_alloc_recfac(int l_max);

/*
 subroutine do_lam_lm(lmax, m, cth, sth, mfac, recfac, lam_lm)
    !=======================================================================
    ! computes scalar lambda_lm(theta) for all l in [m,lmax] for a given m, and given theta
    ! input: lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        recfac: precomputed (by gen_recfac) quantities useful for 
    !            lambda_lm recursion for a given m
    ! output: lam_lm
    ! the routine also needs the array rescale_tac initialized by init_rescale
    !=======================================================================
    integer(I4B),                    intent(in)  :: lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(    0:lmax), intent(out) :: lam_lm
*/



/*
 !=======================================================================
  subroutine gen_mfac( m_max, m_fact)
  !=======================================================================
    ! generates factor used in lam_mm calculation
    ! for all m in 0<=m<=m_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: m_max
    real(DP),     intent(OUT), dimension(0:m_max)  :: m_fact
*/



/*
  !=======================================================================
  subroutine gen_recfac( l_max, m, recfac)
  !=======================================================================
    ! generates recursion factors used to computes the Ylm of degree m 
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                            :: l_max, m
    real(DP),     intent(OUT), dimension(0:1, 0:l_max)  :: recfac
    !
*/
