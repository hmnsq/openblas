#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
தசம NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   தசம ret;
   PRINT_DEBUG_NAME;

   எனில் (n <= 0) திரும்பு 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   எனில் (incx < 0) x -= (n - 1) * incx;
   எனில் (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   திரும்பு ret;
 }

#இல்லையெனில்

தசம CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  தசம ret;
  PRINT_DEBUG_CNAME;

  எனில் (n <= 0) திரும்பு 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  எனில் (incx < 0) x -= (n - 1) * incx;
  எனில் (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  திரும்பு ret;
}

#endif
