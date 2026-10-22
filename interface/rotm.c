#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#ਜੇ_ਅੰਤ

    ROTM_K(n, dx, incx, dy, incy, dparam);

    ਵਾਪਸ;
}

