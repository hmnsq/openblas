#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    بغ_DEBUG_NAME
#وإلا

فراغ CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    بغ_DEBUG_CNAME;

#endif

    ROTM_K(n, dx, incx, dy, incy, dparam);

    أرجع;
}

