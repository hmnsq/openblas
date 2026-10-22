#شامل <stdio.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
اعشاری NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   اعشاری ret;
   PRINT_DEBUG_NAME;

   اگر (n <= 0) بازگشت 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   اگر (incx < 0) x -= (n - 1) * incx;
   اگر (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   بازگشت ret;
 }

#وگرنه

اعشاری CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  اعشاری ret;
  PRINT_DEBUG_CNAME;

  اگر (n <= 0) بازگشت 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  اگر (incx < 0) x -= (n - 1) * incx;
  اگر (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  بازگشت ret;
}

#endif
