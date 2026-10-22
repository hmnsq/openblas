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

   ᬬᬾᬦ᭄ (n <= 0) ᬯᬮᬶ 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   ᬬᬾᬦ᭄ (incx < 0) x -= (n - 1) * incx;
   ᬬᬾᬦ᭄ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   ᬯᬮᬶ ret;
 }

#ᬮᬶᬬ

float CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  float ret;
  PRINT_DEBUG_CNAME;

  ᬬᬾᬦ᭄ (n <= 0) ᬯᬮᬶ 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ᬬᬾᬦ᭄ (incx < 0) x -= (n - 1) * incx;
  ᬬᬾᬦ᭄ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  ᬯᬮᬶ ret;
}

#endif
