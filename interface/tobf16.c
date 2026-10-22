#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#만약 defined(DOUBLE_PREC)
#define FLOAT_TYPE double
#elif defined(SINGLE_PREC)
#define FLOAT_TYPE float
#아니면
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   만약 (n <= 0) 반환;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   만약 (inc_in < 0)   in -= (n - 1) * inc_in;
   만약 (inc_out < 0) out -= (n - 1) * inc_out;

#만약 defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#아니면
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#아니면
void CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  만약 (n <= 0) 반환;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  만약 (inc_in < 0)   in -= (n - 1) * inc_in;
  만약 (inc_out < 0) out -= (n - 1) * inc_out;

#만약 defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
