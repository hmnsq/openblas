#ضمن "common.h"
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#وإلا

فراغ CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#جيڪڏهن_پورو

    ROTM_K(n, dx, incx, dy, incy, dparam);

    أرجع;
}

