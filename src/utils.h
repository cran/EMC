/*
 *  $Id: utils.h,v 1.10 2008/01/11 02:31:17 goswami Exp $
 *  
 *  File:    utils.h
 *  Package: EMC
 *
 *  Created by Gopi Goswami on Wed Apr 12 2006
 *  Copyright (C) 2006 Gopika R. Goswami
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  For a copy of the GNU General Public License please write to the
 *  Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330.
 *  Boston, MA  02111-1307 USA.
 *
 *  For bugs in the code please contact:
 *  goswami@stat.harvard.edu
 *
 *
 *  SYNOPSIS
 *
 *  This file is the declaration file for utils.c; it also contains
 *  some useful macros.
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include <R.h>
#include <Rinternals.h>

/*
 * The #defines
 */
#define MAX_WORD_LENGTH 1024
#define MAX_LINE_LENGTH 4096

#define SQR(_xx) ((_xx) * (_xx))
#define MIN(_aa, _bb) ((_aa) < (_bb) ? (_aa) : (_bb));
#define MAX(_aa, _bb) ((_aa) > (_bb) ? (_aa) : (_bb));
#define PRINT_STUB_INT(_var) Rprintf(#_var " = %d\n", _var);
#define PRINT_STUB_DOUBLE(_var) Rprintf(#_var " = %g\n", _var);
#define PRINT_STUB_STRING(_var) Rprintf(#_var " = %s\n", _var);

#define DEBUG(_xx) \
do { \
      Rprintf("===== BEGIN: DEBUG block =====\n" \
               "file:  %s,  line:  %d\n", \
               __FILE__, __LINE__); \
       _xx \
       Rprintf("===== E N D: DEBUG block =====\n"); \
} while (0)

#define PHONY(_xx) ((void) 0)

#define SWAP(_type, _obj1, _obj2) \
do { \
_type _objTmp; \
\
_objTmp = _obj2; \
_obj2 = _obj1; \
_obj1 = _objTmp; \
} while (0)


extern int
utils_iarray_print (int *arr, int nn, char *sep);

extern int
utils_darray_print (double *arr, int nn, char *sep);

extern int
utils_sarray_print (char **arr, int nn, char *sep);

extern int
utils_SEXP_darray_print (SEXP arr, char *sep);

#endif /* UTILS_H */

