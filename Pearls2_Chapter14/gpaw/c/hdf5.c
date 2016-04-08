/*
 *  Copyright (C) 2010-2011     CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

/* Light weight Python interface to HDF5 functions needed by GPAW

   Generally, HDF5 object identifiers are integers and they are
   passed as such between Python and C. All other data to 
   HDF5 functions is passed in NumPy arrays. At the moment no checks are 
   made for ensuring that data actually is NumPy array.

   The Python interface functions return either
     object identifier as int
     None
   with the exceptions
     tuple from h5s_get_shape
     string from h5l_get_name_by_idx 
*/

   

#ifdef GPAW_WITH_HDF5

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include <hdf5.h>
#include "extensions.h"
#ifdef PARALLEL
#include <mpi.h>
#include "mympi.h"
#endif

// File functions
PyObject* h5f_open(PyObject *self, PyObject *args)
{
  int pid = H5P_DEFAULT;
  const char* name;
  const char mode = 'r';
  if (!PyArg_ParseTuple(args, "s|ci", &name, &mode, &pid))
    return NULL;

  unsigned flag;
  if (mode == 'r')
    flag = H5F_ACC_RDONLY;
  else
    flag = H5F_ACC_RDWR;

  int fid = H5Fopen(name, flag, pid);
  return Py_BuildValue("i", fid);
}

PyObject* h5f_create(PyObject *self, PyObject *args)
{
  int pid = H5P_DEFAULT;  // Property list id
  const char* name;
  if (!PyArg_ParseTuple(args, "s|i", &name, &pid))
    return NULL;

  unsigned flag = H5F_ACC_TRUNC; // Always truncate the file
  int fid = H5Fcreate(name, flag, H5P_DEFAULT, pid);
  return Py_BuildValue("i", fid);
}

PyObject* h5f_close(PyObject *self, PyObject *args)
{
  int fid;
  if (!PyArg_ParseTuple(args, "i", &fid))
    return NULL;

  H5Fclose(fid);
  Py_RETURN_NONE;
}

// Group functions
PyObject* h5g_open(PyObject *self, PyObject *args)
{
  int loc_id;  
  const char* name;
  if (!PyArg_ParseTuple(args, "is", &loc_id, &name))
    return NULL;

  int gid = H5Gopen(loc_id, name, H5P_DEFAULT);
  return Py_BuildValue("i", gid);
}

PyObject* h5g_create(PyObject *self, PyObject *args)
{
  int loc_id;  
  const char* name;
  if (!PyArg_ParseTuple(args, "is", &loc_id, &name))
    return NULL;

  int gid = H5Gcreate2(loc_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  return Py_BuildValue("i", gid);
}

PyObject* h5g_get_num_objs(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5G_info_t group_info;
  H5Gget_info(id, &group_info);
  int nobjs = group_info.nlinks;
  return Py_BuildValue("i", nobjs);
}
  


PyObject* h5g_close(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5Gclose(id);
  Py_RETURN_NONE;
}

// Attribute functions

PyObject* h5a_open(PyObject *self, PyObject *args)
{
  int loc_id;  
  const char* name;
  if (!PyArg_ParseTuple(args, "is", &loc_id, &name))
    return NULL;

  // Check first for the existence
  if (H5Aexists(loc_id, name) == 0 )
    return PyErr_Format(PyExc_KeyError, "HDF5 Attribute %s does not  exist",
                        name);
  int aid = H5Aopen(loc_id, name, H5P_DEFAULT);
  return Py_BuildValue("i", aid);
}

PyObject* h5a_create(PyObject *self, PyObject *args)
{
  int loc_id;  
  int datatype;
  int dataspace;
  const char* name;
  if (!PyArg_ParseTuple(args, "isii", &loc_id, &name, &datatype, &dataspace))
    return NULL;

  hid_t aid = H5Acreate2(loc_id, name, datatype, dataspace, 
		         H5P_DEFAULT, H5P_DEFAULT);

  return Py_BuildValue("i", aid);
}

PyObject* h5a_write(PyObject *self, PyObject *args)
{
  int aid;  
  int datatype;
  PyArrayObject* data;
  if (!PyArg_ParseTuple(args, "iiO", &aid, &datatype, &data))
    return NULL;

  char* buf = PyArray_DATA(data);

  H5Awrite(aid, datatype, buf); 

  Py_RETURN_NONE;
}

PyObject* h5a_read(PyObject *self, PyObject *args)
{
  int aid;  
  int datatype;
  PyArrayObject* data;
  if (!PyArg_ParseTuple(args, "iiO", &aid, &datatype, &data))
    return NULL;

  char* buf = PyArray_DATA(data);

  H5Aread(aid, datatype, buf); 

  Py_RETURN_NONE;
}

PyObject* h5a_get_space(PyObject *self, PyObject *args)
{
  // Returns the dataspace as tuple
  int aid;  
  if (!PyArg_ParseTuple(args, "i", &aid))
    return NULL;

  hid_t dataspace = H5Aget_space(aid); 
  return Py_BuildValue("i", dataspace);
}

PyObject* h5a_get_type(PyObject *self, PyObject *args)
{
  // Returns the dataspace as tuple
  int aid;  
  if (!PyArg_ParseTuple(args, "i", &aid))
    return NULL;

  hid_t filetype = H5Aget_type(aid); 
  hid_t datatype = H5Tget_native_type(filetype, H5T_DIR_ASCEND);
  H5Tclose(filetype);
  return Py_BuildValue("i", datatype);
}

PyObject* h5a_exists_by_name(PyObject *self, PyObject *args)
{
  int id;
  const char* name;
  if (!PyArg_ParseTuple(args, "is", &id, &name))
    return NULL;

  int exists = H5Aexists_by_name(id, ".", name, H5P_DEFAULT);
  return Py_BuildValue("b", exists);
}

PyObject* h5a_close(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5Aclose(id);
  Py_RETURN_NONE;
}

// Datatype functions
PyObject* h5t_get_class(PyObject *self, PyObject *args)
{
  int tid;  
  if (!PyArg_ParseTuple(args, "i", &tid))
    return NULL;

  hid_t class = H5Tget_class(tid); 
  return Py_BuildValue("i", class);
}

PyObject* h5t_get_size(PyObject *self, PyObject *args)
{
  int tid;  
  if (!PyArg_ParseTuple(args, "i", &tid))
    return NULL;

  hid_t size = H5Tget_size(tid); 
  return Py_BuildValue("i", size);
}

PyObject* h5_type_from_numpy(PyObject *self, PyObject *args)
{
  PyArrayObject *array;
  if (!PyArg_ParseTuple(args, "O", &array))
    return NULL;

  int type = PyArray_TYPE(array);
  hid_t datatype;
  if (type == NPY_STRING ) {
    datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, PyArray_ITEMSIZE(array));
  } else if (type == NPY_DOUBLE) {
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  } else if (type == NPY_LONG) {
    datatype = H5Tcopy(H5T_NATIVE_LONG);
  } else if (type == NPY_INT) {
    datatype = H5Tcopy(H5T_NATIVE_INT);
  } else if (type == NPY_BOOL) {
    datatype = H5Tenum_create(H5T_NATIVE_INT8);
    int value;
    value = 0;
    // Convert the int value to int8 with HDF5
    H5Tconvert(H5T_NATIVE_INT, H5T_NATIVE_INT8, 1, &value, NULL, H5P_DEFAULT);
    H5Tenum_insert(datatype, "FALSE", &value);
    value = 1;
    H5Tconvert(H5T_NATIVE_INT, H5T_NATIVE_INT8, 1, &value, NULL, H5P_DEFAULT);
    H5Tenum_insert(datatype, "TRUE", &value);
  } else if (type == NPY_CDOUBLE) {
    datatype = H5Tcreate(H5T_COMPOUND, sizeof(double complex));
    H5Tinsert(datatype, "re", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype, "im", sizeof(double), H5T_NATIVE_DOUBLE);
  } else {
    return PyErr_Format(PyExc_RuntimeError, "Unsupportted datatype");
  }

  return Py_BuildValue("i", datatype);
}

PyObject* h5t_close(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5Tclose(id);
  Py_RETURN_NONE;
}

// Dataspace functions
PyObject* h5s_create(PyObject *self, PyObject *args)
{

  PyArrayObject *shape;
  if (!PyArg_ParseTuple(args, "O", &shape))
    return NULL;

  int rank = PyArray_DIM(shape, 0);
  long* dims_i = PyArray_DATA(shape);
  // hsize_t may be larger than long so we need to copy
  hsize_t* dims = (hsize_t *) malloc(rank * sizeof(hsize_t));
  for (int i=0; i < rank; i++)
    dims[i] = dims_i[i];
   

  int sid = H5Screate_simple(rank, dims, NULL);
  free(dims);
  return Py_BuildValue("i", sid);
}

PyObject* h5s_select_hyperslab(PyObject *self, PyObject *args)
{
  int dataspace;
  PyArrayObject* np_offset;
  PyArrayObject* np_stride;
  PyArrayObject* np_count;
  PyArrayObject* np_block;
  if (!PyArg_ParseTuple(args, "iOOOO", &dataspace, &np_offset, &np_stride, 
			&np_count, &np_block))
    return NULL;
  
  // None can be passed to indicate use of default values e.g. NULL for
  // stride and block
  long* temp = PyArray_DATA(np_offset);
  int rank = PyArray_DIMS(np_offset)[0];
  hsize_t* offset = (hsize_t *) malloc(rank * sizeof(hsize_t));
  for (int i=0; i < rank; i++)
    offset[i] = temp[i];
  hsize_t* stride = NULL;
  if ((PyObject *)np_stride != Py_None) 
    {
      temp = PyArray_DATA(np_stride);
      stride = (hsize_t *) malloc(rank * sizeof(hsize_t));
      for (int i=0; i < rank; i++)
        stride[i] = temp[i];
    }
  temp = PyArray_DATA(np_count); 
  hsize_t* count = (hsize_t *) malloc(rank * sizeof(hsize_t));
  for (int i=0; i < rank; i++)
    count[i] = temp[i]; 
  hsize_t* block = NULL;
  if ((PyObject *)np_block != Py_None) 
    return PyErr_Format(PyExc_NotImplementedError, "Block parameter");

  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride, count, block);

  free(offset);
  if (stride != NULL) 
    free(stride);
  free(count);

  Py_RETURN_NONE;
}

PyObject* h5s_select_none(PyObject *self, PyObject *args)
{
  int dataspace;  
  if (!PyArg_ParseTuple(args, "i", &dataspace))
    return NULL;

  H5Sselect_none(dataspace);
  Py_RETURN_NONE;
}

PyObject* h5s_get_shape(PyObject *self, PyObject *args)
{
  // Returns the dataspace as tuple
  int dataspace;  
  if (!PyArg_ParseTuple(args, "i", &dataspace))
    return NULL;

  int rank = H5Sget_simple_extent_ndims(dataspace);
  hsize_t* dims = (hsize_t *) malloc(rank * sizeof(hsize_t));
  rank = H5Sget_simple_extent_dims(dataspace, dims, NULL);
  PyObject* shape = PyTuple_New(rank);
  int i;
  for (i=0; i < rank; i++)
    {
      PyTuple_SetItem(shape, i, PyInt_FromLong(dims[i]));
    }
  free(dims);
  return shape;
}

PyObject* h5s_close(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5Sclose(id);
  Py_RETURN_NONE;
}


// Dataset functions

PyObject* h5d_open(PyObject *self, PyObject *args)
{
  int loc_id;  
  const char* name;
  if (!PyArg_ParseTuple(args, "is", &loc_id, &name))
    return NULL;

  int did = H5Dopen2(loc_id, name, H5P_DEFAULT);
  return Py_BuildValue("i", did);
}


PyObject* h5d_create(PyObject *self, PyObject *args)
{
  int loc_id;  
  const char* name;
  int dtype_id;
  int space_id;
  if (!PyArg_ParseTuple(args, "isii", &loc_id, &name, &dtype_id, &space_id))
    return NULL;

  int did = H5Dcreate2(loc_id, name, dtype_id, space_id, 
		      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  return Py_BuildValue("i", did);
}

PyObject* h5d_write(PyObject *self, PyObject *args)
{
  int did;  
  int memtype;
  int memspace;  
  int filespace;
  int pid = H5P_DEFAULT;
  PyArrayObject* data;
  if (!PyArg_ParseTuple(args, "iiiiO|i", &did, &memtype, &memspace, 
			&filespace, &data, &pid))
    return NULL;

  char* buf = PyArray_DATA(data);
  H5Dwrite(did, memtype, memspace, filespace, pid, buf);
  
  Py_RETURN_NONE;
}

PyObject* h5d_read(PyObject *self, PyObject *args)
{
  int did;  
  int memtype;
  int memspace;  
  int filespace;
  int pid = H5P_DEFAULT;
  PyArrayObject* data;
  if (!PyArg_ParseTuple(args, "iiiiO|i", &did, &memtype, &memspace, 
			&filespace, &data, &pid))
    return NULL;

  char* buf = PyArray_DATA(data);

  H5Dread(did, memtype, memspace, filespace, pid, buf);
  
  Py_RETURN_NONE;
}

PyObject* h5d_get_space(PyObject *self, PyObject *args)
{
  int did;  
  if (!PyArg_ParseTuple(args, "i", &did))
    return NULL;

  hid_t dataspace = H5Dget_space(did); 
  return Py_BuildValue("i", dataspace);
}

PyObject* h5d_get_type(PyObject *self, PyObject *args)
{
  int did;  
  if (!PyArg_ParseTuple(args, "i", &did))
    return NULL;

  hid_t filetype = H5Dget_type(did); 
  hid_t datatype = H5Tget_native_type(filetype, H5T_DIR_ASCEND);
  H5Tclose(filetype);
  return Py_BuildValue("i", datatype);
}

PyObject* h5d_close(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5Dclose(id);
  Py_RETURN_NONE;
}

// Property list related functions
PyObject* h5p_create(PyObject *self, PyObject *args)
{
  int cls_id;  
  if (!PyArg_ParseTuple(args, "i", &cls_id))
    return NULL;

  int pid = H5Pcreate(cls_id);
  return Py_BuildValue("i", pid);
}

#ifdef PARALLEL
PyObject* h5p_set_fapl_mpio(PyObject *self, PyObject *args)
{
  PyObject *comm_obj;
  int plist_id;
  if (!PyArg_ParseTuple(args, "iO", &plist_id, &comm_obj))
    return NULL;

  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Info info = MPI_INFO_NULL;
  if (comm_obj != Py_None)
    {
      comm = ((MPIObject*)comm_obj)->comm;
      // The following was needed at some point due to bug in Cray MPI
      /* int nprocs;   
      MPI_Comm_size(comm, &nprocs);
      char tmp[20];
      MPI_Info_create(&info);
      sprintf(tmp,"%d", nprocs);
      MPI_Info_set(info,"cb_nodes",tmp); */
#ifdef __bgp__
      // Wavefunction write requires large amounts of memory.                   
      // Appears to be a deficiency in the ROMIO driver.
      // bgl_nodes_pset controls the number of aggregator
      // tasks per pset. The default value is 8.
      // In many cases, it will also be necessary to set
      // DCMF_ALLTOALL_PREMALLOC = N
      char tmp[20];
      int info_param = 32;
      MPI_Info_create(&info);
      sprintf(tmp,"%d", info_param);
      MPI_Info_set(info,"bgl_nodes_pset",tmp);
#endif
    }
  H5Pset_fapl_mpio(plist_id, comm, info);
  Py_RETURN_NONE;
}

PyObject* h5p_set_dxpl_mpio(PyObject *self, PyObject *args)
{
  int plist_id;
  if (!PyArg_ParseTuple(args, "i", &plist_id))
    return NULL;

  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  Py_RETURN_NONE;
}
#endif

PyObject* h5p_close(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;

  H5Pclose(id);
  Py_RETURN_NONE;
}


// Info functions
PyObject* h5i_get_type(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;
  
  int type = H5Iget_type(id);
  return Py_BuildValue("i", type);
}

// Object functions
PyObject* h5o_open(PyObject *self, PyObject *args)
{
  int loc_id;
  const char* name;
  if (!PyArg_ParseTuple(args, "is", &loc_id, &name))
    return NULL;

  int oid = H5Oopen(loc_id, name, H5P_DEFAULT);
  return Py_BuildValue("i", oid);
}

PyObject* h5o_close(PyObject *self, PyObject *args)
{
  int oid;
  if (!PyArg_ParseTuple(args, "i", &oid))
    return NULL;

  H5Oclose(oid);
  Py_RETURN_NONE;
}

// List functions
PyObject* h5l_get_name_by_idx(PyObject *self, PyObject *args)
{
  int id;
  int idx;
  if (!PyArg_ParseTuple(args, "ii", &id, &idx))
    return NULL;

  ssize_t size;
  char *name;
  
  // Get size of the name, add 1 for NULL terminator
  size = 1 + H5Lget_name_by_idx(id, ".", H5_INDEX_NAME, H5_ITER_INC,
				idx, NULL, 0, H5P_DEFAULT);
  name = (char *) malloc(size);
  H5Lget_name_by_idx(id, ".", H5_INDEX_NAME, H5_ITER_INC, idx, name,
		     (size_t) size, H5P_DEFAULT);

  PyObject* retval = Py_BuildValue("s", name);
  free(name);
  return retval;
}


static PyMethodDef functions[] = {
  {"h5f_open", h5f_open, METH_VARARGS, 0},
  {"h5f_create", h5f_create, METH_VARARGS, 0},
  {"h5f_close", h5f_close, METH_VARARGS, 0},
  {"h5g_open", h5g_open, METH_VARARGS, 0},
  {"h5g_get_num_objs", h5g_get_num_objs, METH_VARARGS, 0},
  {"h5g_create", h5g_create, METH_VARARGS, 0},
  {"h5g_close", h5g_close, METH_VARARGS, 0},
  {"h5a_open", h5a_open, METH_VARARGS, 0},
  {"h5a_create", h5a_create, METH_VARARGS, 0},
  {"h5a_write", h5a_write, METH_VARARGS, 0},
  {"h5a_read", h5a_read, METH_VARARGS, 0},
  {"h5a_get_space", h5a_get_space, METH_VARARGS, 0},
  {"h5a_get_type", h5a_get_type, METH_VARARGS, 0},
  {"h5a_exists_by_name", h5a_exists_by_name, METH_VARARGS, 0},
  {"h5a_close", h5a_close, METH_VARARGS, 0},
  {"h5_type_from_numpy", h5_type_from_numpy, METH_VARARGS, 0},
  {"h5t_get_class", h5t_get_class, METH_VARARGS, 0},
  {"h5t_get_size", h5t_get_size, METH_VARARGS, 0},
  {"h5t_close", h5t_close, METH_VARARGS, 0},
  {"h5s_create", h5s_create, METH_VARARGS, 0},
  {"h5s_select_hyperslab", h5s_select_hyperslab, METH_VARARGS, 0},
  {"h5s_select_none", h5s_select_none, METH_VARARGS, 0},
  {"h5s_get_shape", h5s_get_shape, METH_VARARGS, 0},
  {"h5s_close", h5s_close, METH_VARARGS, 0},
  {"h5d_open", h5d_open, METH_VARARGS, 0},
  {"h5d_create", h5d_create, METH_VARARGS, 0},
  {"h5d_write", h5d_write, METH_VARARGS, 0},
  {"h5d_read", h5d_read, METH_VARARGS, 0},
  {"h5d_get_space", h5d_get_space, METH_VARARGS, 0},
  {"h5d_get_type", h5d_get_type, METH_VARARGS, 0},
  {"h5d_close", h5d_close, METH_VARARGS, 0}, 
  {"h5p_create", h5p_create, METH_VARARGS, 0}, 
#ifdef PARALLEL
  {"h5p_set_fapl_mpio", h5p_set_fapl_mpio, METH_VARARGS, 0}, 
  {"h5p_set_dxpl_mpio", h5p_set_dxpl_mpio, METH_VARARGS, 0}, 
#endif
  {"h5p_close", h5p_close, METH_VARARGS, 0}, 
  {"h5o_open", h5o_open, METH_VARARGS, 0}, 
  {"h5o_close", h5o_close, METH_VARARGS, 0}, 
  {"h5i_get_type", h5i_get_type, METH_VARARGS, 0}, 
  {"h5l_get_name_by_idx", h5l_get_name_by_idx, METH_VARARGS, 0},
  {0, 0, 0, 0}
};

PyMODINIT_FUNC init_gpaw_hdf5(void) 
{ 
  PyObject *m = Py_InitModule("_gpaw_hdf5",functions); 
  // Set some hdf5 constants as attributes
  PyModule_AddIntConstant(m, "H5T_FLOAT", H5T_FLOAT);
  PyModule_AddIntConstant(m, "H5T_INTEGER", H5T_INTEGER);
  PyModule_AddIntConstant(m, "H5T_COMPOUND", H5T_COMPOUND);
  PyModule_AddIntConstant(m, "H5T_STRING", H5T_STRING);
  PyModule_AddIntConstant(m, "H5T_ENUM", H5T_ENUM);
  PyModule_AddIntConstant(m, "H5P_DATASET_XFER", H5P_DATASET_XFER);
  PyModule_AddIntConstant(m, "H5P_FILE_ACCESS", H5P_FILE_ACCESS);
  PyModule_AddIntConstant(m, "H5P_DEFAULT", H5P_DEFAULT);
  PyModule_AddIntConstant(m, "H5I_GROUP", H5I_GROUP);
  PyModule_AddIntConstant(m, "H5I_DATASET", H5I_DATASET);
} 
#endif
