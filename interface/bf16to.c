#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ꦪꦼꦤ꧀ defined(DOUBLE_PREC)
#define FLOAT_TYPE dviguu_na
#elif defined(SINGLE_PREC)
#define FLOAT_TYPE float
#ꦭꦶꦪ
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
void NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  ꦪꦼꦤ꧀ (n <= 0) ꦮꦭꦶ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ꦪꦼꦤ꧀ (inc_in < 0)   in -= (n - 1) * inc_in;
  ꦪꦼꦤ꧀ (inc_out < 0) out -= (n - 1) * inc_out;

#ꦪꦼꦤ꧀ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#ꦭꦶꦪ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#ꦭꦶꦪ
void CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  ꦪꦼꦤ꧀ (n <= 0) ꦮꦭꦶ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ꦪꦼꦤ꧀ (inc_in < 0)   in -= (n - 1) * inc_in;
  ꦪꦼꦤ꧀ (inc_out < 0) out -= (n - 1) * inc_out;

#ꦪꦼꦤ꧀ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#ꦭꦶꦪ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
