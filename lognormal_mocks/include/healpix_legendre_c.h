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


/////////////////////////////////////////
//
// C Wrappers for Healpix Fortran Library Assoc. Legendre Functions
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


/*  See the healpix code for detailed documentation on the meaning of the variables in:
      subroutine do_lam_lm(lmax, m, cth, sth, mfac, recfac, lam_lm)
      subroutine gen_mfac( m_max, m_fact)
      subroutine gen_recfac( l_max, m, recfac)
*/
