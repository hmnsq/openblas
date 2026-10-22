#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

OPENBLAS_EXPORT
#اگر_غیر_تعریف_شدہ CBLAS
اعشاری NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   اعشاری ret;
   PRINT_DEBUG_NAME;

   اگر (n <= 0) واپس 0.;

   IDEBUG_آغاز;
   FUNCTION_PROFILE_آغاز();

   اگر (incx < 0) x -= (n - 1) * incx;
   اگر (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_انتہا(1, 2 * n, 2 * n);
   IDEBUG_انتہا;

   واپس ret;
 }

#ورنہ

اعشاری CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  اعشاری ret;
  PRINT_DEBUG_CNAME;

  اگر (n <= 0) واپس 0.;

  IDEBUG_آغاز;
  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0) x -= (n - 1) * incx;
  اگر (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_انتہا(1, 2 * n, 2 * n);
  IDEBUG_انتہا;

  واپس ret;
}

#اگر_ختم
