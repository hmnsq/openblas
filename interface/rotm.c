#ضمن "common.h"
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#وإلا

فراغ CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#د_که_پای

    ROTM_K(n, dx, incx, dy, incy, dparam);

    أرجع;
}

