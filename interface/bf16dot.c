#ضمن <stdio.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
عائم NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   عائم ret;
   PRINT_DEBUG_NAME;

   אם (n <= 0) החזר 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   אם (incx < 0) x -= (n - 1) * incx;
   אם (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   החזר ret;
 }

#אחרת

عائم CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  عائم ret;
  PRINT_DEBUG_CNAME;

  אם (n <= 0) החזר 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  אם (incx < 0) x -= (n - 1) * incx;
  אם (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  החזר ret;
}

#endif
