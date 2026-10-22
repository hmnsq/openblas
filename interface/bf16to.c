#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ಯದಿ defined(DOUBLE_PREC)
#define FLOAT_TYPE ದ್ವಿಗುಣ
#ನೋಚೇದ್ಯದಿ defined(SINGLE_PREC)
#define FLOAT_TYPE ಭಿನ್ನ
#ಅನ್ಯಥಾ
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
ಶೂನ್ಯ NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ಯದಿ (inc_in < 0)   in -= (n - 1) * inc_in;
  ಯದಿ (inc_out < 0) out -= (n - 1) * inc_out;

#ಯದಿ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#ನೋಚೇದ್ಯದಿ defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#ಅನ್ಯಥಾ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#ಅನ್ಯಥಾ
ಶೂನ್ಯ CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ಯದಿ (inc_in < 0)   in -= (n - 1) * inc_in;
  ಯದಿ (inc_out < 0) out -= (n - 1) * inc_out;

#ಯದಿ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#ನೋಚೇದ್ಯದಿ defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#ಅನ್ಯಥಾ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
