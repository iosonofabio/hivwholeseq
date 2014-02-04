/*
 *  Created on:	26/01/2014
 *  Author:	Fabio Zanini
 *  Contents:	Test implementation file
 */
%module seqanpy

%{
#define SWIG_FILE_WITH_INIT
#include "seqanpy.h"
%}

/* STL types */
%include "std_string.i"

/* SEQANPY (move to a separate file?) */
%typemap(in, numinputs=0) int *scoreout(int temp) {
    $1 = &temp;
}
%typemap(argout) int *scoreout {
    PyObject *alipy = PyInt_FromSize_t(*($1));
    PyObject *pos_seed = $result;

    $result = PyTuple_New(2);
    PyTuple_SetItem($result, 0, pos_seed);
    PyTuple_SetItem($result, 1, alipy); 

}


%include "seqanpy.h"
