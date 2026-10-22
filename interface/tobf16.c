#شامل <stdio.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#اگر defined(DOUBLE_PREC)
#define FLOAT_TYPE دوگانه
#وگرنهاگر defined(SINGLE_PREC)
#define FLOAT_TYPE اعشاری
#وگرنه
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
تهی NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   اگر (n <= 0) بازگشت;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   اگر (inc_in < 0)   in -= (n - 1) * inc_in;
   اگر (inc_out < 0) out -= (n - 1) * inc_out;

#اگر defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#وگرنهاگر defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#وگرنه
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#وگرنه
تهی CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  اگر (n <= 0) بازگشت;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  اگر (inc_in < 0)   in -= (n - 1) * inc_in;
  اگر (inc_out < 0) out -= (n - 1) * inc_out;

#اگر defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#وگرنهاگر defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
