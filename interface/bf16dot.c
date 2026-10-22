#ضمن <stdio.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

OPENBLAS_EXPORT
#اگر_تعریف_نشده CBLAS
عائم NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   عائم ret;
   PRINT_DEBUG_NAME;

   إذا (n <= 0) أرجع 0.;

   IDEBUG_شروع;
   FUNCTION_PROFILE_شروع();

   إذا (incx < 0) x -= (n - 1) * incx;
   إذا (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_انتها(1, 2 * n, 2 * n);
   IDEBUG_انتها;

   أرجع ret;
 }

#وإلا

عائم CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  عائم ret;
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0.;

  IDEBUG_شروع;
  FUNCTION_PROFILE_شروع();

  إذا (incx < 0) x -= (n - 1) * incx;
  إذا (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_انتها(1, 2 * n, 2 * n);
  IDEBUG_انتها;

  أرجع ret;
}

#پایان_اگر
