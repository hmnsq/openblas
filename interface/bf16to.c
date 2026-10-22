#ضمن <stdio.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ܐܢ defined(DOUBLE_PREC)
#define FLOAT_TYPE مزدوج
#وإلاإذا defined(SINGLE_PREC)
#define FLOAT_TYPE عائم
#ܐܚܪܢܐ
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
فراغ NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  ܐܢ (n <= 0) ܦܢܐ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ܐܢ (inc_in < 0)   in -= (n - 1) * inc_in;
  ܐܢ (inc_out < 0) out -= (n - 1) * inc_out;

#ܐܢ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#وإلاإذا defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#ܐܚܪܢܐ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#ܐܚܪܢܐ
فراغ CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  ܐܢ (n <= 0) ܦܢܐ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ܐܢ (inc_in < 0)   in -= (n - 1) * inc_in;
  ܐܢ (inc_out < 0) out -= (n - 1) * inc_out;

#ܐܢ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#وإلاإذا defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#ܐܚܪܢܐ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
