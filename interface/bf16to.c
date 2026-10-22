#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ប្រសិនបើ defined(DOUBLE_PREC)
#define FLOAT_TYPE dviguu_na
#elif defined(SINGLE_PREC)
#define FLOAT_TYPE float
#បើមិនដូច្នេះ
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
void NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  ប្រសិនបើ (n <= 0) ត្រឡប់;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ប្រសិនបើ (inc_in < 0)   in -= (n - 1) * inc_in;
  ប្រសិនបើ (inc_out < 0) out -= (n - 1) * inc_out;

#ប្រសិនបើ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#បើមិនដូច្នេះ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#បើមិនដូច្នេះ
void CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  ប្រសិនបើ (n <= 0) ត្រឡប់;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ប្រសិនបើ (inc_in < 0)   in -= (n - 1) * inc_in;
  ប្រសិនបើ (inc_out < 0) out -= (n - 1) * inc_out;

#ប្រសិនបើ defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#elif defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#បើមិនដូច្នេះ
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
