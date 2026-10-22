#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#ບໍ່ດັ່ງນັ້ນ

void CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#endif

    ROTM_K(n, dx, incx, dy, incy, dparam);

    ສົ່ງຄືນ;
}

