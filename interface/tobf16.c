#ضمن <stdio.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#𐡄𐡍 defined(DOUBLE_PREC)
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

   PRINT_DEBUG_NAME;

   𐡄𐡍 (n <= 0) أرجع;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   𐡄𐡍 (inc_in < 0)   in -= (n - 1) * inc_in;
   𐡄𐡍 (inc_out < 0) out -= (n - 1) * inc_out;

#𐡄𐡍 defined(DOUBLE_PREC)
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
  PRINT_DEBUG_CNAME;

  𐡄𐡍 (n <= 0) أرجع;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  𐡄𐡍 (inc_in < 0)   in -= (n - 1) * inc_in;
  𐡄𐡍 (inc_out < 0) out -= (n - 1) * inc_out;

#𐡄𐡍 defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#وإلاإذا defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
