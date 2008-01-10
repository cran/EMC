/*
 *  $Id: utils.c,v 1.4 2008/01/11 02:31:16 goswami Exp $
 *  
 *  File:    utils.c
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
 *  This file contains all the utility functions used by the functions
 *  in the file objects.c.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include "utils.h"

int
utils_iarray_print (int *arr, int nn, char *sep)
{
	int ii;

	if ((arr == NULL) || (nn < 0))
	  Rprintf("MALFORMED array\n");
	else if (nn == 0)
	  Rprintf("EMPTY array\n");
	else {
		for (ii = 0; ii < (nn - 1); ++ii) 
		  Rprintf("%d%s", arr[ii], sep);
		Rprintf("%d\n", arr[ii]);
	}
        return 0;
}


int
utils_darray_print (double *arr, int nn, char *sep)
{
	int ii;

	if ((arr == NULL) || (nn < 0))
	  Rprintf("MALFORMED array\n");
	else if (nn == 0)
	  Rprintf("EMPTY array\n");
	else {
		for (ii = 0; ii < (nn - 1); ++ii) 
		  Rprintf("%g%s", arr[ii], sep);
		Rprintf("%g\n", arr[ii]);
	}
        return 0;
}


int
utils_sarray_print (char **arr, int nn, char *sep)
{
	int ii;

	if ((arr == NULL) || (nn < 0))
	  Rprintf("MALFORMED array\n");
	else if (nn == 0)
	  Rprintf("EMPTY array\n");
	else {
		for (ii = 0; ii < (nn - 1); ++ii) 
		  Rprintf("%s%s", arr[ii], sep);
		Rprintf("%s\n", arr[ii]);
	}
        return 0;
}


int
utils_SEXP_darray_print (SEXP arr, char *sep)
{
	int ii, nn = length(arr);

	if ((arr == NULL) || (nn < 0))
	  Rprintf("MALFORMED array\n");
	else if (nn == 0)
	  Rprintf("EMPTY array\n");
	else {
		for (ii = 0; ii < (nn - 1); ++ii) 
		  Rprintf("%g%s", REAL(arr)[ii], sep);
		Rprintf("%g\n", REAL(arr)[ii]);
	}
        return 0;
}


