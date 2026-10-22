#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
float NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   float ret;
   PRINT_DEBUG_NAME;

   ถ้า (n <= 0) คืน 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   ถ้า (incx < 0) x -= (n - 1) * incx;
   ถ้า (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   คืน ret;
 }

#มิฉะนั้น

float CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  float ret;
  PRINT_DEBUG_CNAME;

  ถ้า (n <= 0) คืน 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ถ้า (incx < 0) x -= (n - 1) * incx;
  ถ้า (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  คืน ret;
}

#endif
