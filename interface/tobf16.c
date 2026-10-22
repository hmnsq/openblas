#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#यदि defined(DOUBLE_PREC)
#define FLOAT_TYPE द्विगुण
#अथवायदि defined(SINGLE_PREC)
#define FLOAT_TYPE भिन्न
#अन्यथा
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
शून्य NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   यदि (n <= 0) फर्काऊ;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   यदि (inc_in < 0)   in -= (n - 1) * inc_in;
   यदि (inc_out < 0) out -= (n - 1) * inc_out;

#यदि defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#अथवायदि defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#अन्यथा
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#अन्यथा
शून्य CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  यदि (n <= 0) फर्काऊ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  यदि (inc_in < 0)   in -= (n - 1) * inc_in;
  यदि (inc_out < 0) out -= (n - 1) * inc_out;

#यदि defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#अथवायदि defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
