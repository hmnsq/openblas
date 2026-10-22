#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
ভিন্ন NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   ভিন্ন ret;
   PRINT_DEBUG_NAME;

   যদি (n <= 0) প্রত্যাগম 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   যদি (incx < 0) x -= (n - 1) * incx;
   যদি (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   প্রত্যাগম ret;
 }

#অন্যথা

ভিন্ন CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  ভিন্ন ret;
  PRINT_DEBUG_CNAME;

  যদি (n <= 0) প্রত্যাগম 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (n - 1) * incx;
  যদি (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  প্রত্যাগম ret;
}

#endif
