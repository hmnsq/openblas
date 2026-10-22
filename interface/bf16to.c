#ضمن <stdio.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#إذا defined(DOUBLE_PREC)
#تعریف_کن FLOAT_TYPE مزدوج
#وإلاإذا defined(SINGLE_PREC)
#تعریف_کن FLOAT_TYPE عائم
#وإلا
#پایان_اگر

OPENBLAS_EXPORT
#اگر_تعریف_نشده CBLAS
فراغ NAME(blasint *N, bfloat16 *درون, blasint *INC_درون, FLOAT_TYPE *بیرون, blasint *INC_بیرون){
  BLASLONG n    = *N;
  BLASLONG inc_درون = *INC_درون;
  BLASLONG inc_بیرون = *INC_بیرون;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع;

  IDEBUG_شروع;
  FUNCTION_PROFILE_شروع();

  إذا (inc_درون < 0)   درون -= (n - 1) * inc_درون;
  إذا (inc_بیرون < 0) بیرون -= (n - 1) * inc_بیرون;

#إذا defined(DOUBLE_PREC)
  D_BF16_مقصد_K(n, درون, inc_درون, بیرون, inc_بیرون);
#وإلاإذا defined(SINGLE_PREC)
  S_BF16_مقصد_K(n, درون, inc_درون, بیرون, inc_بیرون);
#وإلا
#پایان_اگر

  FUNCTION_PROFILE_انتها(1, 2 * n, 2 * n);
  IDEBUG_انتها;
}
#وإلا
فراغ CNAME(blasint n, bfloat16 * درون, blasint inc_درون, FLOAT_TYPE * بیرون, blasint inc_بیرون){
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع;

  IDEBUG_شروع;
  FUNCTION_PROFILE_شروع();

  إذا (inc_درون < 0)   درون -= (n - 1) * inc_درون;
  إذا (inc_بیرون < 0) بیرون -= (n - 1) * inc_بیرون;

#إذا defined(DOUBLE_PREC)
  D_BF16_مقصد_K(n, درون, inc_درون, بیرون, inc_بیرون);
#وإلاإذا defined(SINGLE_PREC)
  S_BF16_مقصد_K(n, درون, inc_درون, بیرون, inc_بیرون);
#وإلا
#پایان_اگر

  FUNCTION_PROFILE_انتها(1, 2 * n, 2 * n);
  IDEBUG_انتها;
}
#پایان_اگر
