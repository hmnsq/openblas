#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
દશાંશ NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   દશાંશ ret;
   PRINT_DEBUG_NAME;

   જો (n <= 0) પાછા 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   જો (incx < 0) x -= (n - 1) * incx;
   જો (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   પાછા ret;
 }

#નહીંતર

દશાંશ CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  દશાંશ ret;
  PRINT_DEBUG_CNAME;

  જો (n <= 0) પાછા 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  જો (incx < 0) x -= (n - 1) * incx;
  જો (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  પાછા ret;
}

#endif
