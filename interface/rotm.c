#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(blasint *N, FLOAT *dx, blasint *INCX, FLOAT *dy, blasint *INCY, FLOAT *dparam){

    blasint n = *N;
    blasint incx = *INCX;
    blasint incy = *INCY;

    PRINT_DEBUG_NAME
#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(blasint n, FLOAT *dx, blasint incx, FLOAT *dy, blasint incy, FLOAT *dparam){

    PRINT_DEBUG_CNAME;

#ಆದರೆ_ಅಂತ್ಯ

    ROTM_K(n, dx, incx, dy, incy, dparam);

    ಹಿಂತಿರುಗಿ;
}

