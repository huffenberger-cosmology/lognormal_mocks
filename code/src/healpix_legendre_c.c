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




#include <healpix_legendre_c.h>

double *hpl_alloc_recfac(int l_max){
 
  return( (double *) calloc(2*(l_max+1),sizeof(double)) );
}

double *hpl_alloc_mfac(int m_max){

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
