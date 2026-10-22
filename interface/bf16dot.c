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

   सचे (n <= 0) पटिनिवत्त 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   सचे (incx < 0) x -= (n - 1) * incx;
   सचे (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   पटिनिवत्त ret;
 }

#अञ्ञथा

भिन्न CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  भिन्न ret;
  PRINT_DEBUG_CNAME;

  सचे (n <= 0) पटिनिवत्त 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  सचे (incx < 0) x -= (n - 1) * incx;
  सचे (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  पटिनिवत्त ret;
}

#endif
