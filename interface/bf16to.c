#ضمن <stdio.h>
#ضمن "common.h"
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#إذا defined(DOUBLE_PREC)
#تعریف_کړه FLOAT_TYPE مزدوج
#وإلاإذا defined(SINGLE_PREC)
#تعریف_کړه FLOAT_TYPE عائم
#وإلا
#د_که_پای

OPENBLAS_EXPORT
#که_نه_تعریف_شوی CBLAS
فراغ NAME(blasint *N, bfloat16 *دننه, blasint *INC_دننه, FLOAT_TYPE *بهر, blasint *INC_بهر){
  BLASLONG n    = *N;
  BLASLONG inc_دننه = *INC_دننه;
  BLASLONG inc_بهر = *INC_بهر;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع;

  IDEBUG_پیل;
  FUNCTION_PROFILE_پیل();

  إذا (inc_دننه < 0)   دننه -= (n - 1) * inc_دننه;
  إذا (inc_بهر < 0) بهر -= (n - 1) * inc_بهر;

#إذا defined(DOUBLE_PREC)
  D_BF16_منزل_K(n, دننه, inc_دننه, بهر, inc_بهر);
#وإلاإذا defined(SINGLE_PREC)
  S_BF16_منزل_K(n, دننه, inc_دننه, بهر, inc_بهر);
#وإلا
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(1, 2 * n, 2 * n);
  IDEBUG_پای_ټکی;
}
#وإلا
فراغ CNAME(blasint n, bfloat16 * دننه, blasint inc_دننه, FLOAT_TYPE * بهر, blasint inc_بهر){
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع;

  IDEBUG_پیل;
  FUNCTION_PROFILE_پیل();

  إذا (inc_دننه < 0)   دننه -= (n - 1) * inc_دننه;
  إذا (inc_بهر < 0) بهر -= (n - 1) * inc_بهر;

#إذا defined(DOUBLE_PREC)
  D_BF16_منزل_K(n, دننه, inc_دننه, بهر, inc_بهر);
#وإلاإذا defined(SINGLE_PREC)
  S_BF16_منزل_K(n, دننه, inc_دننه, بهر, inc_بهر);
#وإلا
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(1, 2 * n, 2 * n);
  IDEBUG_پای_ټکی;
}
#د_که_پای
