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
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>


void W2Cl(double *theta, double *W, int Ntheta, double *l, double *C, int Nl);
void Cl2W(double *theta, double *W, int Ntheta, double *l, double *C, int Nl);
void Cl2W_nospline(double *theta, double *W, int Ntheta, double *C, int Nl); // uses every ell
