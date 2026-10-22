#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
ദശാംശം NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   ദശാംശം ret;
   PRINT_DEBUG_NAME;

   എങ്കിൽ (n <= 0) തിരികെ 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   എങ്കിൽ (incx < 0) x -= (n - 1) * incx;
   എങ്കിൽ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   തിരികെ ret;
 }

#അല്ലെങ്കിൽ

ദശാംശം CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  ദശാംശം ret;
  PRINT_DEBUG_CNAME;

  എങ്കിൽ (n <= 0) തിരികെ 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0) x -= (n - 1) * incx;
  എങ്കിൽ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  തിരികെ ret;
}

#endif
