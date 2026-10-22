#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
ਦਸ਼ਮਲਵ NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   ਦਸ਼ਮਲਵ ret;
   PRINT_DEBUG_NAME;

   ਜੇ (n <= 0) ਵਾਪਸ 0.;

   IDEBUG_ਸ਼ੁਰੂ;
   FUNCTION_PROFILE_ਸ਼ੁਰੂ();

   ਜੇ (incx < 0) x -= (n - 1) * incx;
   ਜੇ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_ਸਿਰਾ(1, 2 * n, 2 * n);
   IDEBUG_ਸਿਰਾ;

   ਵਾਪਸ ret;
 }

#ਨਹੀਂਤਾਂ

ਦਸ਼ਮਲਵ CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  ਦਸ਼ਮਲਵ ret;
  PRINT_DEBUG_CNAME;

  ਜੇ (n <= 0) ਵਾਪਸ 0.;

  IDEBUG_ਸ਼ੁਰੂ;
  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0) x -= (n - 1) * incx;
  ਜੇ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_ਸਿਰਾ(1, 2 * n, 2 * n);
  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;
}

#ਜੇ_ਅੰਤ
