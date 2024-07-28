#include <float.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <arrayobject.h>

#include "sym_assert.h"
#include "types.h"
#include "gen/snavely_reprojection_factor.h"
#include "gen/snavely_reprojection_factor2.h"
#include "gen/rot3_tangent.h"

static PyObject *
read_problem(PyObject *self, PyObject *args) {
   char* path = NULL;

    if (!PyArg_ParseTuple(args, 
                            "s", 
                           &path)) {
        return NULL;
    }

    f64 epsilon = 10.0 * DBL_EPSILON;

    FILE* file = fopen(path, "r");

    i32 num_cameras, num_points, num_observations;
    fscanf(file, "%d", &num_cameras);
    fscanf(file, "%d", &num_points);
    fscanf(file, "%d", &num_observations);

    npy_intp camera_indices_dim[1] = {num_observations};
    PyObject* camera_indices = PyArray_New(&PyArray_Type, 1, camera_indices_dim, NPY_INT32, NULL, NULL, 0, NPY_ARRAY_IN_ARRAY, NULL);
    i32* camera_indices_data = (i32*) PyArray_DATA(camera_indices);

    npy_intp point_indices_dim[1] = {num_observations};
    PyObject* point_indices = PyArray_New(&PyArray_Type, 1, point_indices_dim, NPY_INT32, NULL, NULL, 0, NPY_ARRAY_IN_ARRAY, NULL);
    i32* point_indices_data = (i32*) PyArray_DATA(point_indices);

    npy_intp pixels_dim[2] = {num_observations, 2};
    PyObject* pixels = PyArray_New(&PyArray_Type, 2, pixels_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_ARRAY, NULL);
    f64* pixels_data = (f64*) PyArray_DATA(pixels);

    for (i32 i = 0; i < num_observations; i++) {
        i32 camera, point;
        fscanf(file, "%d", &camera);
        fscanf(file, "%d", &point);

        f64 px, py;
        fscanf(file, "%lf", &px);
        fscanf(file, "%lf", &py);

        camera_indices_data[i] = camera;
        point_indices_data[i] = point;
        pixels_data[2 * i + 0] = px;
        pixels_data[2 * i + 1] = py;
    }

    // For some reason NPY_ARRAY_IN_ARRAY doesn't give a C-contiguous array, but 0 does.

    // (cam_T_world: 7, intrinsics: 3, point: 3)
    npy_intp cam_T_world_dim[2] = {num_cameras, 7};
    PyObject* cam_T_world = PyArray_New(&PyArray_Type, 2, cam_T_world_dim, NPY_DOUBLE, NULL, NULL, 0, 0, NULL);
    SYM_ASSERT(PyArray_ISALIGNED(cam_T_world));
    SYM_ASSERT(PyArray_IS_C_CONTIGUOUS(cam_T_world));
    f64* cam_T_world_data = (f64*) PyArray_DATA(cam_T_world);

    npy_intp intrinsics_dim[2] = {num_cameras, 3};
    PyObject* intrinsics = PyArray_New(&PyArray_Type, 2, intrinsics_dim, NPY_DOUBLE, NULL, NULL, 0, 0, NULL);
    SYM_ASSERT(PyArray_ISALIGNED(intrinsics));
    SYM_ASSERT(PyArray_IS_C_CONTIGUOUS(intrinsics));
    f64* intrinsics_data = (f64*) PyArray_BYTES(intrinsics);

    for (i32 i = 0; i < num_cameras; i++) {
        f64 rx, ry, rz, tx, ty, tz, f, k1, k2;
        fscanf(file, "%lf", &rx);
        fscanf(file, "%lf", &ry);
        fscanf(file, "%lf", &rz);
        fscanf(file, "%lf", &tx);
        fscanf(file, "%lf", &ty);
        fscanf(file, "%lf", &tz);
        fscanf(file, "%lf", &f);
        fscanf(file, "%lf", &k1);
        fscanf(file, "%lf", &k2);

        f64 rot_tangent[3] = {rx, ry, rz};
        rot3_from_tangent(rot_tangent, &cam_T_world_data[7 * i], epsilon);

        cam_T_world_data[7 * i + 4] = tx;
        cam_T_world_data[7 * i + 5] = ty;
        cam_T_world_data[7 * i + 6] = tz;

        intrinsics_data[3 * i + 0] = f;
        intrinsics_data[3 * i + 1] = k1;
        intrinsics_data[3 * i + 2] = k2;
    }

    npy_intp point_dim[2] = {num_points, 3};
    PyObject* point = PyArray_New(&PyArray_Type, 2, point_dim, NPY_DOUBLE, NULL, NULL, 0, 0, NULL);
    SYM_ASSERT(PyArray_ISALIGNED(point));
    SYM_ASSERT(PyArray_IS_C_CONTIGUOUS(point));
    f64* point_data = (f64*) PyArray_DATA(point);

    for (i32 i = 0; i < num_points; i++) {
        f64 x, y, z;
        fscanf(file, "%lf", &x);
        fscanf(file, "%lf", &y);
        fscanf(file, "%lf", &z);

        point_data[3 * i + 0] = x;
        point_data[3 * i + 1] = y;
        point_data[3 * i + 2] = z;
    }

    fclose(file);

    return Py_BuildValue("NNNNNN", cam_T_world, intrinsics, point, camera_indices, point_indices, pixels);
}

#define SYM_DOUBLE_VEC(m) \
        PyArray_FromAny(m, PyArray_DescrFromType(NPY_DOUBLE), 1, 1, NPY_ARRAY_IN_FARRAY, NULL)

static PyObject *
snavely_reprojection_factor_py(PyObject *self, PyObject *args)
{
    PyObject* cam_T_world_storage_arg = NULL;
    PyObject* intrinsics_storage_arg = NULL;
    PyObject* point_storage_arg = NULL;
    PyObject* pixel_storage_arg = NULL;

    f64 epsilon = 10.0 * DBL_EPSILON;

    if (!PyArg_ParseTuple(args, 
                            "OOOO", 
                            &cam_T_world_storage_arg,
                            &intrinsics_storage_arg, 
                            &point_storage_arg, 
                            &pixel_storage_arg)) {
        return NULL;
    }

    PyObject* cam_T_world_storage = SYM_DOUBLE_VEC(cam_T_world_storage_arg);
    SYM_ASSERT(PyArray_DIM(cam_T_world_storage, 0) == 7);

    PyObject* intrinsics_storage = SYM_DOUBLE_VEC(intrinsics_storage_arg);
    SYM_ASSERT(PyArray_DIM(intrinsics_storage, 0) == 3);

    PyObject* point_storage = SYM_DOUBLE_VEC(point_storage_arg);
    SYM_ASSERT(PyArray_DIM(point_storage, 0) == 3);

    PyObject* pixel_storage = SYM_DOUBLE_VEC(pixel_storage_arg);
    SYM_ASSERT(PyArray_DIM(pixel_storage, 0) == 2);


    npy_intp res_dim[1] = {2};
    PyObject* res_storage = PyArray_New(&PyArray_Type, 1, res_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_FARRAY, NULL);
    PyArray_FILLWBYTE(res_storage, 0);

    npy_intp hessian_dim[2] = {12, 12};
    PyObject* hessian_storage = PyArray_New(&PyArray_Type, 2, hessian_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_FARRAY, NULL);

    npy_intp rhs_dim[1] = {12};
    PyObject* rhs_storage = PyArray_New(&PyArray_Type, 1, rhs_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_FARRAY, NULL);

    // Jacobians are broken!
    snavely_reprojection_factor(PyArray_DATA(cam_T_world_storage),
                                 PyArray_DATA(intrinsics_storage),
                                 PyArray_DATA(point_storage),
                                 PyArray_DATA(pixel_storage),
                                 epsilon, PyArray_DATA(res_storage), 
                                    /*jacobian_storage*/ NULL, 
                                    PyArray_DATA(hessian_storage),
                                PyArray_DATA(rhs_storage));

    return Py_BuildValue("NNN", res_storage, hessian_storage, rhs_storage);
}

static PyObject *
snavely_reprojection_factor2_py(PyObject *self, PyObject *args)
{
    PyObject* cam_T_world_storage_arg = NULL;
    PyObject* intrinsics_storage_arg = NULL;
    PyObject* point_storage_arg = NULL;
    PyObject* pixel_storage_arg = NULL;

    f64 epsilon = 10.0 * DBL_EPSILON;

    if (!PyArg_ParseTuple(args, 
                            "OOOO", 
                            &cam_T_world_storage_arg,
                            &intrinsics_storage_arg, 
                            &point_storage_arg, 
                            &pixel_storage_arg)) {
        return NULL;
    }

    PyObject* cam_T_world_storage = SYM_DOUBLE_VEC(cam_T_world_storage_arg);
    SYM_ASSERT(PyArray_DIM(cam_T_world_storage, 0) == 7);

    PyObject* intrinsics_storage = SYM_DOUBLE_VEC(intrinsics_storage_arg);
    SYM_ASSERT(PyArray_DIM(intrinsics_storage, 0) == 3);

    PyObject* point_storage = SYM_DOUBLE_VEC(point_storage_arg);
    SYM_ASSERT(PyArray_DIM(point_storage, 0) == 3);

    PyObject* pixel_storage = SYM_DOUBLE_VEC(pixel_storage_arg);
    SYM_ASSERT(PyArray_DIM(pixel_storage, 0) == 2);


    npy_intp res_dim[1] = {2};
    PyObject* res_storage = PyArray_New(&PyArray_Type, 1, res_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_FARRAY, NULL);
    PyArray_FILLWBYTE(res_storage, 0);

    npy_intp hessian_dim[2] = {12, 12};
    PyObject* hessian_storage = PyArray_New(&PyArray_Type, 2, hessian_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_FARRAY, NULL);

    npy_intp rhs_dim[1] = {12};
    PyObject* rhs_storage = PyArray_New(&PyArray_Type, 1, rhs_dim, NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_IN_FARRAY, NULL);

    // Jacobians are broken!
    snavely_reprojection_factor2(PyArray_DATA(cam_T_world_storage),
                                 PyArray_DATA(intrinsics_storage),
                                 PyArray_DATA(point_storage),
                                 PyArray_DATA(pixel_storage),
                                 epsilon, PyArray_DATA(res_storage), 
                                    /*jacobian_storage*/ NULL, 
                                    PyArray_DATA(hessian_storage),
                                PyArray_DATA(rhs_storage));

    return Py_BuildValue("NNN", res_storage, hessian_storage, rhs_storage);
}

static PyMethodDef BalMethods[] = {
    {"snavely_reprojection_factor",  snavely_reprojection_factor_py, METH_VARARGS,
     "Compute a Snavely reprojection factor."},
    {"snavely_reprojection_factor2",  snavely_reprojection_factor2_py, METH_VARARGS,
     "Compute a Snavely reprojection factor."},
    {"read_problem",  read_problem, METH_VARARGS,
     "Read a BAL problem from a file."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef balmodule = {
    PyModuleDef_HEAD_INIT,
    "bal",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    BalMethods
};

PyMODINIT_FUNC
PyInit_bal(void)
{
    import_array(); // PyError if not successful
    return PyModule_Create(&balmodule);
}

