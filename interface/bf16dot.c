#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
ଦଶମିକ NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   ଦଶମିକ ret;
   PRINT_DEBUG_NAME;

   ଯଦି (n <= 0) ଫେରନ୍ତୁ 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   ଯଦି (incx < 0) x -= (n - 1) * incx;
   ଯଦି (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   ଫେରନ୍ତୁ ret;
 }

#ଅନ୍ୟଥା

ଦଶମିକ CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  ଦଶମିକ ret;
  PRINT_DEBUG_CNAME;

  ଯଦି (n <= 0) ଫେରନ୍ତୁ 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0) x -= (n - 1) * incx;
  ଯଦି (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  ଫେରନ୍ତୁ ret;
}

#endif
