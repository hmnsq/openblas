#ضمن <stdio.h>
#ضمن "common.h"
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

OPENBLAS_EXPORT
#که_نه_تعریف_شوی CBLAS
عائم NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   عائم ret;
   PRINT_DEBUG_NAME;

   إذا (n <= 0) أرجع 0.;

   IDEBUG_پیل;
   FUNCTION_PROFILE_پیل();

   إذا (incx < 0) x -= (n - 1) * incx;
   إذا (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_پای_ټکی(1, 2 * n, 2 * n);
   IDEBUG_پای_ټکی;

   أرجع ret;
 }

#وإلا

عائم CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  عائم ret;
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0.;

  IDEBUG_پیل;
  FUNCTION_PROFILE_پیل();

  إذا (incx < 0) x -= (n - 1) * incx;
  إذا (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_پای_ټکی(1, 2 * n, 2 * n);
  IDEBUG_پای_ټکی;

  أرجع ret;
}

#د_که_پای
