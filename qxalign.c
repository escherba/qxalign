/*
 * =====================================================================================
 *
 *       Filename:  qxalign.c
 *
 *    Description:  Python 3 module exposing quality-aware alignment routines in
 *                  align454.c, align454.h
 *
 *        Version:  1.0
 *        Created:  05/07/2011 17:45:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eugene Scherba (es), escherba@bu.edu
 *        Company:  Boston University
 *
 * =====================================================================================
 */

#include <Python.h>
#include "structmember.h"
#include "align454.h"

/*-----------------------------------------------------------------------------
 *  Qxalign type object
 *-----------------------------------------------------------------------------*/
typedef struct {
        PyObject_HEAD
        Alignment_ASW* al;
        Py_buffer db_seq;
        Py_buffer query_seq;
        Py_buffer query_qual;
        uint8_t* default_qual;
        int match;
        int mismatch;
        int gap_open_extend;
        int gap_extend;
} Qxalign;

/*-----------------------------------------------------------------------------
 *  Custom exception objects
 *-----------------------------------------------------------------------------*/
/* static PyObject *QxalignError; */

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_new
 *  Description:  tp_new
 * =====================================================================================
 */
static PyObject *
Qxalign_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
        Qxalign *self;

        if ((self = (Qxalign *)type->tp_alloc(type, 0)) == NULL) {
                PyErr_SetString(PyExc_MemoryError, "cannot allocate Qxalign instance");
                return NULL;
        }

        if ((self->al = asw_alloc(PyMem_Malloc, PyMem_Realloc, PyMem_Free)) == NULL) {
                Py_DECREF(self);
                PyErr_SetString(PyExc_MemoryError, "cannot allocate alignment object");
                return NULL;
        }

        /* self->db_seq.buf = NULL; */
        /* self->query_seq.buf = NULL; */
        /* self->query_qual.buf = NULL; */
        self->db_seq.len = 0u;
        self->query_seq.len = 0u;
        self->query_qual.len = 0u;

        /* Set default values for penalty scores */
        self->match = -10;
        self->mismatch = 30;
        self->gap_open_extend = 50;
        self->gap_extend = 20;

        self->default_qual = NULL;

        return (PyObject *)self;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_dealloc
 *  Description:  deallocate an instance of Qxalign
 * =====================================================================================
 */
static void
Qxalign_dealloc(Qxalign* self)
{
        asw_free(self->al);

        PyBuffer_Release(&(self->db_seq));
        PyBuffer_Release(&(self->query_seq));
        PyBuffer_Release(&(self->query_qual));

        if (self->default_qual != NULL) {
                PyMem_Free(self->default_qual);
        }

        Py_TYPE(self)->tp_free((PyObject*)self);
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_init
 *  Description:  tp_init
 * =====================================================================================
 */
static int
Qxalign_init(Qxalign *self, PyObject *args, PyObject *kwds)
{
        static char *kwlist[] =
                {"match", "mismatch", "gap_open_extend", "gap_extend", NULL};

        if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiii", kwlist,
                                &self->match,
                                &self->mismatch,
                                &self->gap_open_extend,
                                &self->gap_extend))
        {
                return -1;
        }
        asw_init(self->al,
                 self->match,
                 self->mismatch,
                 self->gap_open_extend,
                 self->gap_extend);
        return 0;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_prepare_db
 *  Description:  Parse arguments and resize alignment object
 * =====================================================================================
 */
static PyObject*
Qxalign_prepare_db(Qxalign* self, PyObject *args, PyObject *kwds)
{
        Py_buffer db_seq;
        db_seq.buf = NULL;

        static char *kwlist[] = {
                "db_seq",
                NULL /*  Sentinel */
        };
        if (!PyArg_ParseTupleAndKeywords(args, kwds, "z*", kwlist,
                                         &db_seq))
        {
                return NULL;
        }
        if (db_seq.buf != NULL) {
                PyBuffer_Release(&(self->db_seq));
                self->db_seq = db_seq;
        }
        if (asw_prepare_db(self->al,
                    (const char*)self->db_seq.buf,
                    self->db_seq.len, 0u, 0u)
                != 0)
        {
                PyErr_SetString(PyExc_MemoryError, "cannot resize alignment object");
                return NULL;
        }
        Py_RETURN_NONE;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_prepare_query
 *  Description:  Parse arguments and resize alignment object
 * =====================================================================================
 */
static PyObject*
Qxalign_prepare_query(Qxalign* self, PyObject *args, PyObject *kwds)
{
        Py_buffer query_seq, query_qual;
        query_seq.buf = NULL;
        query_qual.buf = NULL;

        int phred_offset = 33,
            assume_phred = PHRED_RANGE - 1;

        static char *kwlist[] = {
                "query_seq",
                "query_qual",
                "phred_offset",
                "assume_phred",
                NULL /*  Sentinel */
        };
        if (!PyArg_ParseTupleAndKeywords(args, kwds, "z*|z*ii", kwlist,
                                         &query_seq,
                                         &query_qual,
                                         &phred_offset,
                                         &assume_phred))
        {
                return NULL;
        }
        if (query_seq.buf != NULL) {
                PyBuffer_Release(&(self->query_seq));
                self->query_seq = query_seq;
        }

        uint8_t *final_qual;
        if (query_qual.buf != NULL) {
                if (query_qual.len < self->query_seq.len) {
                        PyErr_SetString(PyExc_IndexError,
                                "quality score array is shorter than query sequence");
                        return NULL;
                } else if (query_qual.len > self->query_seq.len &&
                           PyErr_WarnEx(PyExc_UserWarning,
                            "quality score array is longer than query sequence", 2) < 0)
                {
                        return NULL;
                }
                PyBuffer_Release(&(self->query_qual));
                self->query_qual = query_qual;
                final_qual = query_qual.buf;
        } else {
                if (assume_phred >= PHRED_RANGE || assume_phred < 0) {
                        PyErr_Format(PyExc_IndexError,
                                "assumed PHRED score %d is outside of valid range 0-%d",
                                assume_phred, PHRED_RANGE - 1);
                        return NULL;
                }
                /* maintain a buffer with default PHRED scores */
                uint8_t* tmp = PyMem_Realloc(self->default_qual,
                                             sizeof(uint8_t) * self->query_seq.len);
                if (tmp == NULL) {
                        PyErr_SetString(PyExc_MemoryError,
                                "cannot resize default quality score array");
                        return NULL;
                }
                self->default_qual = tmp;

                memset(tmp, assume_phred + phred_offset, self->query_seq.len);
                final_qual = tmp;
        }

        if (asw_prepare(self->al,
                    (const char*)self->db_seq.buf,
                    self->db_seq.len,
                    (const char*)self->query_seq.buf,
                    (const uint8_t*)final_qual,
                    self->query_seq.len, 0u, 0u)
                != 0)
        {
                PyErr_SetString(PyExc_MemoryError, "cannot resize alignment object");
                return NULL;
        }
        asw_set_phoffset(self->al, phred_offset);

        Py_RETURN_NONE;
}
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_prepare
 *  Description:  Parse arguments and resize alignment object
 * =====================================================================================
 */
static PyObject*
Qxalign_prepare(Qxalign* self, PyObject *args, PyObject *kwds)
{
        Py_buffer db_seq, query_seq, query_qual;

        db_seq.buf = NULL;
        query_seq.buf = NULL;
        query_qual.buf = NULL;

        int phred_offset = 33,
            assume_phred = PHRED_RANGE - 1;

        static char *kwlist[] = {
                "db_seq",
                "query_seq",
                "query_qual",
                "phred_offset",
                "assume_phred",
                NULL /*  Sentinel */
        };
        if (!PyArg_ParseTupleAndKeywords(args, kwds, "z*z*|z*ii", kwlist,
                                         &db_seq,
                                         &query_seq,
                                         &query_qual,
                                         &phred_offset,
                                         &assume_phred))
        {
                return NULL;
        }
        if (db_seq.buf != NULL) {
                PyBuffer_Release(&(self->db_seq));
                self->db_seq = db_seq;
        }
        if (query_seq.buf != NULL) {
                PyBuffer_Release(&(self->query_seq));
                self->query_seq = query_seq;
        }

        uint8_t *final_qual;
        if (query_qual.buf != NULL) {
                if (query_qual.len < self->query_seq.len) {
                        PyErr_SetString(PyExc_IndexError,
                                "quality score array is shorter than query sequence");
                        return NULL;
                } else if (query_qual.len > self->query_seq.len &&
                           PyErr_WarnEx(PyExc_UserWarning,
                            "quality score array is longer than query sequence", 2) < 0)
                {
                        return NULL;
                }
                PyBuffer_Release(&(self->query_qual));
                self->query_qual = query_qual;
                final_qual = query_qual.buf;
        } else {
                if (assume_phred >= PHRED_RANGE || assume_phred < 0) {
                        PyErr_Format(PyExc_IndexError,
                                "assumed PHRED score %d is outside of valid range 0-%d",
                                assume_phred, PHRED_RANGE - 1);
                        return NULL;
                }
                /* maintain a buffer with default PHRED scores */
                uint8_t* tmp = PyMem_Realloc(self->default_qual,
                                             sizeof(uint8_t) * self->query_seq.len);
                if (tmp == NULL) {
                        PyErr_SetString(PyExc_MemoryError,
                                "cannot resize default quality score array");
                        return NULL;
                }
                self->default_qual = tmp;

                memset(tmp, assume_phred + phred_offset, self->query_seq.len);
                final_qual = tmp;
        }

        if (asw_prepare(self->al,
                    (const char*)self->db_seq.buf,
                    self->db_seq.len,
                    (const char*)self->query_seq.buf,
                    (const uint8_t*)final_qual,
                    self->query_seq.len, 0u, 0u)
                != 0)
        {
                PyErr_SetString(PyExc_MemoryError, "cannot resize alignment object");
                return NULL;
        }
        asw_set_phoffset(self->al, phred_offset);

        Py_RETURN_NONE;
}
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_align
 *  Description:  Perofrm alignment
 * =====================================================================================
 */
static PyObject *
Qxalign_align(Qxalign* self, PyObject *args, PyObject *kwds)
{
        PyObject *x = Py_False;
        static char *kwlist[] = {"semi", NULL};
        if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O:bool", kwlist, &x)) {
                return NULL;
        }
        int semi = PyObject_IsTrue(x);

        if (self->al->subdb_len == 0u || self->al->subquery_len == 0u) {
                PyErr_SetString(PyExc_IndexError,
                        "cannot perform alignment on a zero-element matrix");
                return NULL;
        }
        if (semi) {
                /* semiglobal alignment: fill top row (parallel to db) with zeroes */
                asw_align_init_semi(self->al);
        } else {
                /* global alignment: penalize deletions in db at the beginning */
                asw_align_init(self->al);
        }
        asw_align(self->al);
        /* asw_print_matrix1(self->al, stdout); */
        return Py_BuildValue("i", asw_locate_minscore(self->al));
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_trace
 *  Description:  Performs traceback on an alignment
 * =====================================================================================
 */
static PyObject *
Qxalign_trace(Qxalign* self)
{
        if (self->al->subdb_len == 0u || self->al->subquery_len == 0u) {
                PyErr_SetString(PyExc_IndexError,
                        "cannot perform traceback on a zero-element matrix");
                return NULL;
        }
        asw_trace(self->al);
        Py_RETURN_NONE;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_print_trace
 *  Description:  Print CIGAR traceback of the alignment
 * =====================================================================================
 */
static PyObject *
Qxalign_print_trace(Qxalign* self)
{
        asw_print_cigar(self->al, stdout);
        Py_RETURN_NONE;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Qxalign_show_trace
 *  Description:  Print CIGAR traceback of the alignment
 * =====================================================================================
 */
static PyObject *
Qxalign_show_trace(Qxalign* self)
{
        const char *buf = asw_show_cigar(self->al);
        PyObject *str = PyUnicode_FromString(buf);
        free(buf);
        return str;
}

/*-----------------------------------------------------------------------------
 *  Module-level data fields
 *-----------------------------------------------------------------------------*/
static PyModuleDef qxalign_module = {
        PyModuleDef_HEAD_INIT,
        "qxalign",
        "Quality-aware realignment of sequence reads",
        -1,
        NULL, NULL, NULL, NULL, NULL
};

/*-----------------------------------------------------------------------------
 *  Qxalign type data members
 *-----------------------------------------------------------------------------*/
static PyMemberDef Qxalign_members[] = {
        /* Python string objects */
        /* {"db_seq", T_OBJECT_EX, offsetof(Qxalign, db_seq), 0,
         *        "database (reference) sequence"},
         *{"query_seq", T_OBJECT_EX, offsetof(Qxalign, query_seq), 0,
         *        "query sequence"},
         */

        /* integers */
        {"match", T_INT, offsetof(Qxalign, match), READONLY,
                "match penalty"},
        {"mismatch", T_INT, offsetof(Qxalign, mismatch), READONLY,
                "mismatch penalty"},
        {"gap_open_extend", T_INT, offsetof(Qxalign, gap_open_extend), READONLY,
                "gap open+extend penalty"},
        {"gap_extend", T_INT, offsetof(Qxalign, gap_extend), READONLY,
                "gap extension penalty"},

        {NULL}  /* Sentinel */
};

/*-----------------------------------------------------------------------------
 *  Qxalign type methods
 *-----------------------------------------------------------------------------*/
static PyMethodDef Qxalign_methods[] = {
        {"prepare", (PyCFunction)Qxalign_prepare, METH_VARARGS|METH_KEYWORDS,
                "Assign input sequences and resizes the alignment matrix"},
        {"prepare_db", (PyCFunction)Qxalign_prepare_db, METH_VARARGS|METH_KEYWORDS,
                "Assign reference (database) sequence and resizes the alignment matrix"},
        {"prepare_query", (PyCFunction)Qxalign_prepare_query, METH_VARARGS|METH_KEYWORDS,
                "Assign query sequence, query quality string and resizes the alignment matrix"},
        {"align", (PyCFunction)Qxalign_align, METH_VARARGS|METH_KEYWORDS,
                "Perform an alignment and returns resulting score"},
        {"trace", (PyCFunction)Qxalign_trace, METH_NOARGS,
                "Perform traceback on an alignment"},
        {"print_trace", (PyCFunction)Qxalign_print_trace, METH_NOARGS,
                "Print CIGAR traceback of an alignment to stdout"},
        {"show_trace", (PyCFunction)Qxalign_show_trace, METH_NOARGS,
                "Return CIGAR traceback of an alignment"},
        {NULL}  /* Sentinel */
};

/*-----------------------------------------------------------------------------
 *  Qxalign type type info
 *-----------------------------------------------------------------------------*/
static PyTypeObject QxalignType = {
        PyVarObject_HEAD_INIT(NULL, 0)
                "noddy.Qxalign",             /* tp_name */
        sizeof(Qxalign),             /* tp_basicsize */
        0,                         /* tp_itemsize */
        (destructor)Qxalign_dealloc, /* tp_dealloc */
        0,                         /* tp_print */
        0,                         /* tp_getattr */
        0,                         /* tp_setattr */
        0,                         /* tp_reserved */
        0,                         /* tp_repr */
        0,                         /* tp_as_number */
        0,                         /* tp_as_sequence */
        0,                         /* tp_as_mapping */
        0,                         /* tp_hash  */
        0,                         /* tp_call */
        0,                         /* tp_str */
        0,                         /* tp_getattro */
        0,                         /* tp_setattro */
        0,                         /* tp_as_buffer */
        Py_TPFLAGS_DEFAULT |
                Py_TPFLAGS_BASETYPE,   /* tp_flags */
        "Qxalign objects",           /* tp_doc */
        0,                         /* tp_traverse */
        0,                         /* tp_clear */
        0,                         /* tp_richcompare */
        0,                         /* tp_weaklistoffset */
        0,                         /* tp_iter */
        0,                         /* tp_iternext */
        Qxalign_methods,             /* tp_methods */
        Qxalign_members,             /* tp_members */
        0,                         /* tp_getset */
        0,                         /* tp_base */
        0,                         /* tp_dict */
        0,                         /* tp_descr_get */
        0,                         /* tp_descr_set */
        0,                         /* tp_dictoffset */
        (initproc)Qxalign_init,      /* tp_init */
        0,                         /* tp_alloc */
        Qxalign_new                  /* tp_new */
};

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  PyInit_qxalign
 *  Description:  PyMODINIT_FUNC
 * =====================================================================================
 */
PyMODINIT_FUNC
PyInit_qxalign(void)
{
        PyObject* m;

        /* QxalignType.tp_new = Qxalign_new; */
        if (PyType_Ready(&QxalignType) < 0) {
                return NULL;
        }
        if ((m = PyModule_Create(&qxalign_module)) == NULL) {
                return NULL;
        }

        Py_INCREF(&QxalignType);
        PyModule_AddObject(m, "Qxalign", (PyObject *)&QxalignType);

        /*  create custom exception */
        /* if (QxalignError == NULL) {
         *        QxalignError = PyErr_NewException("qxalign.error", NULL, NULL);
         *}
         *Py_INCREF(QxalignError);
         *PyModule_AddObject(m, "error", QxalignError);
         */

        return m;
}
