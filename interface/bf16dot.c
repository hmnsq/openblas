#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

OPENBLAS_EXPORT
#अगर_अघोषित CBLAS
भग्न NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   भग्न ret;
   PRINT_DEBUG_NAME;

   अगर (n <= 0) वापस 0.;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   अगर (incx < 0) x -= (n - 1) * incx;
   अगर (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;

   वापस ret;
 }

#अन्यथा

भग्न CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  भग्न ret;
  PRINT_DEBUG_CNAME;

  अगर (n <= 0) वापस 0.;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (n - 1) * incx;
  अगर (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;

  वापस ret;
}

#पूर्ण_अगर
