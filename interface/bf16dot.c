#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
ಭಿನ್ನ NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   ಭಿನ್ನ ret;
   PRINT_DEBUG_NAME;

   ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   ಯದಿ (incx < 0) x -= (n - 1) * incx;
   ಯದಿ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   ಪ್ರತ್ಯಾಗಮ ret;
 }

#ಅನ್ಯಥಾ

ಭಿನ್ನ CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  ಭಿನ್ನ ret;
  PRINT_DEBUG_CNAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (n - 1) * incx;
  ಯದಿ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;
}

#endif
