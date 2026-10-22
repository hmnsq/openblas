#ضمن <stdio.h>
#ضمن "common.h"
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

OPENBLAS_EXPORT
#جيڪڏهن_وصف_نه_ٿيل CBLAS
عائم NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   عائم ret;
   PRINT_DEBUG_NAME;

   إذا (n <= 0) أرجع 0.;

   IDEBUG_شروعات;
   FUNCTION_PROFILE_شروعات();

   إذا (incx < 0) x -= (n - 1) * incx;
   إذا (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_پڇاڙي(1, 2 * n, 2 * n);
   IDEBUG_پڇاڙي;

   أرجع ret;
 }

#وإلا

عائم CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  عائم ret;
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0.;

  IDEBUG_شروعات;
  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0) x -= (n - 1) * incx;
  إذا (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_پڇاڙي(1, 2 * n, 2 * n);
  IDEBUG_پڇاڙي;

  أرجع ret;
}

#جيڪڏهن_پورو
