#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#எனில் defined(DOUBLE_PREC)
#define FLOAT_TYPE இரட்டை
#அல்லதெனில் defined(SINGLE_PREC)
#define FLOAT_TYPE தசம
#இல்லையெனில்
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
வெறுமை NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   எனில் (n <= 0) திரும்பு;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   எனில் (inc_in < 0)   in -= (n - 1) * inc_in;
   எனில் (inc_out < 0) out -= (n - 1) * inc_out;

#எனில் defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#அல்லதெனில் defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#இல்லையெனில்
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#இல்லையெனில்
வெறுமை CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  எனில் (n <= 0) திரும்பு;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  எனில் (inc_in < 0)   in -= (n - 1) * inc_in;
  எனில் (inc_out < 0) out -= (n - 1) * inc_out;

#எனில் defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#அல்லதெனில் defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
