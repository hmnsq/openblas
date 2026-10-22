#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

OPENBLAS_EXPORT
#నిర్వచితం_కాకపోతే CBLAS
దశాంశం NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   దశాంశం ret;
   PRINT_DEBUG_NAME;

   ఒకవేళ (n <= 0) తిరిగి 0.;

   IDEBUG_ప్రారంభం;
   FUNCTION_PROFILE_ప్రారంభం();

   ఒకవేళ (incx < 0) x -= (n - 1) * incx;
   ఒకవేళ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_చివర(1, 2 * n, 2 * n);
   IDEBUG_చివర;

   తిరిగి ret;
 }

#లేకపోతే

దశాంశం CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  దశాంశం ret;
  PRINT_DEBUG_CNAME;

  ఒకవేళ (n <= 0) తిరిగి 0.;

  IDEBUG_ప్రారంభం;
  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0) x -= (n - 1) * incx;
  ఒకవేళ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_చివర(1, 2 * n, 2 * n);
  IDEBUG_చివర;

  తిరిగి ret;
}

#అయితే_ముగింపు
