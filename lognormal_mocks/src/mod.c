#include <Python.h>
#include <numpy/ndarrayobject.h>

static PyObject *helloworld(PyObject *self, PyObject *args) {

  printf("mod hello world\n");
  
  return(Py_None);
}


static PyMethodDef lognormal_mocksMethods[] = {
  {"helloworld", helloworld, METH_VARARGS, "hello world\n"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initmod(void)
{
  (void) Py_InitModule("mod", lognormal_mocksMethods);
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
}
