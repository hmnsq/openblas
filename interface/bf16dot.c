#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
भिन्न NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   भिन्न ret;
   PRINT_DEBUG_NAME;

   जइ (n <= 0) प्रत्यागम 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   जइ (incx < 0) x -= (n - 1) * incx;
   जइ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   प्रत्यागम ret;
 }

#अण्णहा

भिन्न CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  भिन्न ret;
  PRINT_DEBUG_CNAME;

  जइ (n <= 0) प्रत्यागम 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  जइ (incx < 0) x -= (n - 1) * incx;
  जइ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  प्रत्यागम ret;
}

#endif
