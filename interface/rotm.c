#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#अन्यथा

शून्य CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#यदि_अन्तः

    ROTM_K(n, dx, incx, dy, incy, dparam);

    प्रत्यागम;
}

