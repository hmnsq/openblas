#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#જો defined(DOUBLE_PREC)
#define FLOAT_TYPE દ્વિગુણ
#અથવાજો defined(SINGLE_PREC)
#define FLOAT_TYPE દશાંશ
#નહીંતર
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
શૂન્ય NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   જો (n <= 0) પાછા;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   જો (inc_in < 0)   in -= (n - 1) * inc_in;
   જો (inc_out < 0) out -= (n - 1) * inc_out;

#જો defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#અથવાજો defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#નહીંતર
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#નહીંતર
શૂન્ય CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  જો (n <= 0) પાછા;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  જો (inc_in < 0)   in -= (n - 1) * inc_in;
  જો (inc_out < 0) out -= (n - 1) * inc_out;

#જો defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#અથવાજો defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
