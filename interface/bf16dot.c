#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
दशांश NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   दशांश ret;
   PRINT_DEBUG_NAME;

   जर (n <= 0) परत 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   जर (incx < 0) x -= (n - 1) * incx;
   जर (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   परत ret;
 }

#अन्यथा

दशांश CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  दशांश ret;
  PRINT_DEBUG_CNAME;

  जर (n <= 0) परत 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  जर (incx < 0) x -= (n - 1) * incx;
  जर (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  परत ret;
}

#endif
