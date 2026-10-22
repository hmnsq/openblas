#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
দশমিক NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   দশমিক ret;
   PRINT_DEBUG_NAME;

   যদি (n <= 0) ফেরত 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   যদি (incx < 0) x -= (n - 1) * incx;
   যদি (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   ফেরত ret;
 }

#অন্যথা

দশমিক CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  দশমিক ret;
  PRINT_DEBUG_CNAME;

  যদি (n <= 0) ফেরত 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (n - 1) * incx;
  যদি (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  ফেরত ret;
}

#endif
