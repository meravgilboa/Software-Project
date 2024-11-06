#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"

/*Decleration of the module methods*/
double** matrix_python_to_c(PyObject *matrix_py, int rows, int cols);
PyObject* matrix_c_to_python(double **matrix_c, int rows, int cols);
static PyObject *symnmf(PyObject *self, PyObject *args);
static PyObject *sym(PyObject *self, PyObject *args);
static PyObject *ddg(PyObject *self, PyObject *args);
static PyObject *norm(PyObject *self, PyObject *args);
static PyObject *derive(PyObject *self, PyObject *args);

double** matrix_python_to_c(PyObject *matrix_py, int rows, int cols){
    PyObject *row;
    double **matrix_c;
    int i, j;
    matrix_c = malloc(rows*sizeof(double*));
    for (i=0; i<rows; i++){
        row = PyList_GetItem(matrix_py, i);
        if (!PyList_Check(row) || PyList_Size(row) != cols)
        {
            PyErr_SetString(PyExc_ValueError, "An Error Has Occured");
            return NULL;
        }
        matrix_c[i] = malloc(cols*sizeof(double));
        for (j=0; j<cols; j++){
            matrix_c[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }
    return matrix_c;
}

PyObject* matrix_c_to_python(double **matrix_c, int rows, int cols){
    PyObject *matrix_py, *row;
    int i, j;
    matrix_py = PyList_New(rows);
    for (i=0; i<rows; i++){
        row = PyList_New(cols);
        for (j=0; j<cols; j++){
            PyList_SetItem(row, j, PyFloat_FromDouble(matrix_c[i][j]));
        }
        PyList_SetItem(matrix_py, i, row);
    }
    return matrix_py;
}

static PyObject *symnmf(PyObject *self, PyObject *args){
    PyObject *W_py, *H_py, *result_py;
    double **W, **H;
    int n, k;
    /*Parse arguments*/
    if (!PyArg_ParseTuple(args, "OOii", &W_py, &H_py, &n, &k))
    {
        return NULL;
    }
    /*Convert to C data structure*/
    W = matrix_python_to_c(W_py, n, n);
    H = matrix_python_to_c(H_py, n, k);
    /*Call the function*/
    double** result = converge_H(H, W, n, k);
    /*Convert to Python data structure*/
    result_py = matrix_c_to_python(result, n, k);
    /*Free memory*/
    free_matrix(W, n);
    /*Return the result*/
    return result_py;
}

static PyObject *sym(PyObject *self, PyObject *args){
    PyObject *datapoints, *sym_mat;
    double **points;
    double **A;
    int n, d;
    /*Parse arguments*/
    if (!PyArg_ParseTuple(args, "Oii", &datapoints, &n, &d))
    {
        return NULL;
    }
    /*Convert to C data structure*/
    points = matrix_python_to_c(datapoints, n, d);
    /*Call the function*/
    A = calc_A(points, n, d);
    /*Convert to Python data structure*/
    sym_mat = matrix_c_to_python(A, n, n);
    /*Free memory*/
    free_matrix(points, n);
    free_matrix(A, n);

    return sym_mat;
}

static PyObject *ddg(PyObject *self, PyObject *args){
    PyObject *datapoints, *diagonal_mat;
    double **points, **A, **D;
    int n, d;
    /*Parse arguments*/
    if (!PyArg_ParseTuple(args, "Oii", &datapoints, &n, &d))
    {
        return NULL;
    }
    /*Convert to C data structure*/
    points = matrix_python_to_c(datapoints, n, d);
    /*Call the function*/
    A = calc_A(points, n, d);
    D = calc_D(A, n);
    /*Convert to Python data structure*/
    diagonal_mat = matrix_c_to_python(D, n, n);
    /*Free memory*/
    free_matrix(points, n);
    free_matrix(A, n);
    free_matrix(D, n);

    return diagonal_mat;
}

static PyObject *norm(PyObject *self, PyObject *args){
    PyObject *datapoints, *norm_mat;
    double **points, **A, **D, **W;
    int n, d;
    /*Parse arguments*/
    if (!PyArg_ParseTuple(args, "Oii", &datapoints, &n, &d))
    {
        return NULL;
    }
    /*Convert to C data structure*/
    points = matrix_python_to_c(datapoints, n, d);
    /*Call the function*/
    A = calc_A(points, n, d);
    D = calc_D(A, n);
    W = calc_W(D, A, n);
    /*Convert to Python data structure*/
    norm_mat = matrix_c_to_python(W, n, n);
    /*Free memory*/
    free_matrix(points, n);
    free_matrix(A, n);
    free_matrix(D, n);

    return norm_mat;
}

static PyObject *derive(PyObject *self, PyObject *args){
    PyObject *H_py, *result_py;
    double **H;
    int n, k, i;
    /*Parse arguments*/
    if (!PyArg_ParseTuple(args, "Oii", &H_py, &n, &k))
    {
        printf("Did not parse arguments\n");
        return NULL;
    }
    /*Convert to C data structure*/
    H = matrix_python_to_c(H_py, n, k);
    /*Call the function*/
    int* result = derive_solution(H, n, k);
    /*Convert to Python data structure*/
    result_py = PyList_New(n);
    for (i=0; i<n; i++)
    {
        PyList_SetItem(result_py, i, PyLong_FromLong(result[i]));
    }
    /*Free memory*/
    free_matrix(H, n);
    free(result);
    
    return result_py;
}
/*Decleration of the module methods*/
static PyMethodDef symnmf_methods[] = {
    {"symnmf", symnmf, METH_VARARGS, PyDoc_STR("Gets W, H, n, k. Executes the symnmf algorithm in C. Returns final H matrix as List of Lists")},
    {"sym", sym, METH_VARARGS, PyDoc_STR("Gets points list, n, dimension. Returns similarity matrix as List of Lists")},
    {"ddg", ddg, METH_VARARGS, PyDoc_STR("Gets points list, n, dimension.Returns diagonal matrix as List of Lists")},
    {"norm", norm, METH_VARARGS, PyDoc_STR("Gets points list, n, dimension.Returns norm matrix as List of Lists")},
    {"derive", derive, METH_VARARGS, PyDoc_STR("Derive the solution of cluster assignments computed in C")},
    {NULL, NULL, 0, NULL}};
/*Module definition*/
static struct PyModuleDef symnmf_module = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symnmf_methods,
};
/*Module initialization function*/
PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *module;
    module = PyModule_Create(&symnmf_module);
    if (module == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    return module;
}