#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#если defined(DOUBLE_PREC)
#define FLOAT_TYPE double
#elif defined(SINGLE_PREC)
#define FLOAT_TYPE float
#иначе
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   если (n <= 0) вернуть;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   если (inc_in < 0)   in -= (n - 1) * inc_in;
   если (inc_out < 0) out -= (n - 1) * inc_out;

#если defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#иначе
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#иначе
void CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  если (n <= 0) вернуть;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  если (inc_in < 0)   in -= (n - 1) * inc_in;
  если (inc_out < 0) out -= (n - 1) * inc_out;

#если defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
