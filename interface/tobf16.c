#ضمن <stdio.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#إذا defined(DOUBLE_PREC)
#define FLOAT_TYPE مزدوج
#وإلاإذا defined(SINGLE_PREC)
#define FLOAT_TYPE عائم
#وإلا
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
فراغ NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   اطبع_DEBUG_NAME;

   إذا (n <= 0) أرجع;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   إذا (inc_in < 0)   in -= (n - 1) * inc_in;
   إذا (inc_out < 0) out -= (n - 1) * inc_out;

#إذا defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#وإلاإذا defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#وإلا
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#وإلا
فراغ CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  اطبع_DEBUG_CNAME;

  إذا (n <= 0) أرجع;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  إذا (inc_in < 0)   in -= (n - 1) * inc_in;
  إذا (inc_out < 0) out -= (n - 1) * inc_out;

#إذا defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#وإلاإذا defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
