#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#якщо defined(DOUBLE_PREC)
#define FLOAT_TYPE double
#elif defined(SINGLE_PREC)
#define FLOAT_TYPE float
#інакше
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
void NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  якщо (n <= 0) повернути;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  якщо (inc_in < 0)   in -= (n - 1) * inc_in;
  якщо (inc_out < 0) out -= (n - 1) * inc_out;

#якщо defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#інакше
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#інакше
void CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  якщо (n <= 0) повернути;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  якщо (inc_in < 0)   in -= (n - 1) * inc_in;
  якщо (inc_out < 0) out -= (n - 1) * inc_out;

#якщо defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#інакше
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
