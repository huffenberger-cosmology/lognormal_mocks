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
  npy_intp *Nlnpy = PyArray_DIMS(Cl_array);
  int Nl = (int) (Nlnpy[2]);

  printf("Nmaps = %d\nNl = %d\n", Nmaps, Nl);
  

  double *gaussbar = calloc(Nmaps,sizeof(double));
  double *Clgauss = calloc(Nmaps*Nmaps*Nl,sizeof(double)); // Power spectrum of Gaussianized field

  
  double *xidelta = calloc(Nmaps*Nmaps*Nth,sizeof(double)); // Correlation function of field

  lognormal_mocks_stats_fullsky(Nmaps,
				rhobar, // mean densities (input)
				Nl,Cl, // power spectrum of field (input)
				gaussbar, Clgauss,  // mean and overdensity spectrum of Gaussianized field (output)
				xidelta, // correlation function of overdensity
				Nth // accuracy parameter
				);

  PyObject *npygaussbar = PyArray_SimpleNewFromData(1, Nmapsnpy, NPY_DOUBLE, gaussbar);
  PyArray_ENABLEFLAGS(npygaussbar, NPY_OWNDATA);

  PyObject *npyClgauss = PyArray_SimpleNewFromData(3, Nlnpy, NPY_DOUBLE, Clgauss);
  PyArray_ENABLEFLAGS(npyClgauss, NPY_OWNDATA);


  PyObject *outtuple = PyTuple_Pack(2, npygaussbar, npyClgauss);
  
  return(outtuple);
}





static PyMethodDef lognormal_mocksMethods[] = {
  {"lognormal_mocks_stats", mod_lognormal_mocks_stats_fullsky, METH_VARARGS, "Returns mean and pverdensity power spectrum of Gaussianized field.  Parameters rhobar[Nmaps], Cl[Nmaps*Nmaps*Nl], Ntheta (accuracy)\n"},
  {"helloworld", helloworld, METH_VARARGS, "hello world\n"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initmod(void)
{
  (void) Py_InitModule("mod", lognormal_mocksMethods);
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
}
