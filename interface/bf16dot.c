#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

OPENBLAS_EXPORT
#यदि_अपरिभाषितम् CBLAS
भिन्न NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   भिन्न ret;
   PRINT_DEBUG_NAME;

   यदि (n <= 0) प्रत्यागम 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   यदि (incx < 0) x -= (n - 1) * incx;
   यदि (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_अन्तः(1, 2 * n, 2 * n);
   IDEBUG_अन्तः;

   प्रत्यागम ret;
 }

#अन्यथा

भिन्न CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  भिन्न ret;
  PRINT_DEBUG_CNAME;

  यदि (n <= 0) प्रत्यागम 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (n - 1) * incx;
  यदि (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_अन्तः(1, 2 * n, 2 * n);
  IDEBUG_अन्तः;

  प्रत्यागम ret;
}

#यदि_अन्तः
