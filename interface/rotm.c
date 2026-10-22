#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#وإلا

فراغ CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#جے_انت

    ROTM_K(n, dx, incx, dy, incy, dparam);

    أرجع;
}

