#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#എങ്കിൽ defined(DOUBLE_PREC)
#define FLOAT_TYPE ഇരട്ട
#അല്ലെങ്കിലെങ്കിൽ defined(SINGLE_PREC)
#define FLOAT_TYPE ദശാംശം
#അല്ലെങ്കിൽ
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
ശൂന്യം NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  എങ്കിൽ (n <= 0) തിരികെ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  എങ്കിൽ (inc_in < 0)   in -= (n - 1) * inc_in;
  എങ്കിൽ (inc_out < 0) out -= (n - 1) * inc_out;

#എങ്കിൽ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#അല്ലെങ്കിലെങ്കിൽ defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#അല്ലെങ്കിൽ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#അല്ലെങ്കിൽ
ശൂന്യം CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  എങ്കിൽ (n <= 0) തിരികെ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  എങ്കിൽ (inc_in < 0)   in -= (n - 1) * inc_in;
  എങ്കിൽ (inc_out < 0) out -= (n - 1) * inc_out;

#എങ്കിൽ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#അല്ലെങ്കിലെങ്കിൽ defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#അല്ലെങ്കിൽ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
