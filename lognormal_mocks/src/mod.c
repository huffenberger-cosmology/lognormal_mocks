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




//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <lognormal_mocks.h>

static PyObject *helloworld(PyObject *self, PyObject *args) {

  printf("mod hello world\n");
  
  return(Py_None);
}



static PyObject *mod_lognormal_mocks_stats_fullsky(PyObject *self, PyObject *args) {

  int Nth;  /* Theta resolution (accuracy parameter) for intermediate correlation function */
  PyObject *rhobar_array;
  PyObject *Cl_array;

  
  if (!PyArg_ParseTuple(args, "OOi", &rhobar_array, &Cl_array,&Nth)) 
    return NULL;

  double *rhobar = PyArray_DATA(rhobar_array);
  double *Cl =  PyArray_DATA(Cl_array);

  npy_intp *Nmapsnpy = PyArray_DIMS(rhobar_array);
  int Nmaps = (int) *PyArray_DIMS(rhobar_array);
  npy_intp *NClnpy = PyArray_DIMS(Cl_array);
  int Nl = (int) (NClnpy[2]);

  printf("Nmaps = %d\nNl = %d\n", Nmaps, Nl);

  if ( (Nmaps != (int)NClnpy[0]) || (Nmaps != (int)NClnpy[1]) ) {
    printf("Problem: incompatible dimensions\n Nmaps = %d NCl = [%d, %d]\n",Nmaps,(int)NClnpy[0],(int)NClnpy[1]);
    return(Py_None);
  }
  
  
  PyArrayObject *npygaussbar = (PyArrayObject *) PyArray_SimpleNew(1, Nmapsnpy,NPY_DOUBLE);
  PyArrayObject *npyClgauss = (PyArrayObject *) PyArray_SimpleNew(3, NClnpy,NPY_DOUBLE);
  
  double *xidelta = calloc(Nmaps*Nmaps*Nth,sizeof(double)); // Correlation function of field

  lognormal_mocks_stats_fullsky(Nmaps,
				rhobar, // mean densities (input)
				Nl,Cl, // power spectrum of field (input)
				(double *)npygaussbar->data, (double *) npyClgauss->data,  // mean and overdensity spectrum of Gaussianized field (output)
				xidelta, // correlation function of overdensity
				Nth // accuracy parameter
				);

  PyObject *tuple = Py_BuildValue("NN", npygaussbar, npyClgauss);

  free(xidelta);
  
  return(tuple);
}





static PyMethodDef lognormal_mocksMethods[] = {
  {"lognormal_mocks_stats", mod_lognormal_mocks_stats_fullsky, METH_VARARGS, "Returns mean and pverdensity power spectrum of Gaussianized field.  Parameters rhobar[Nmaps], Cl[Nmaps*Nmaps*Nl], Ntheta (accuracy)\n"},
  {"helloworld", helloworld, METH_VARARGS, "hello world\n"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

#if PY_MAJOR_VERSION == 2

PyMODINIT_FUNC initmod(void)
{
  (void) Py_InitModule("lognormal_mocks_mod", lognormal_mocksMethods);
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
}
#endif

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef lognormal_mocks_module = {
					     PyModuleDef_HEAD_INIT,
					     "lognormal_mocks_mod",   /* name of module */
					     NULL, /* module documentation, may be NULL */
					     -1,       /* size of per-interpreter state of the module,
							  or -1 if the module keeps state in global variables. */
					     lognormal_mocksMethods
};


PyMODINIT_FUNC PyInit_lognormal_mocks_mod(void)
{
  PyObject *m;
  
  m = PyModule_Create(&lognormal_mocks_module);
  
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
  
  return(m);
}

#endif
