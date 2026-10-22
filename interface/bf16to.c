#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#যদি defined(DOUBLE_PREC)
#define FLOAT_TYPE দ্বিগুণ
#নতুবা_যদি defined(SINGLE_PREC)
#define FLOAT_TYPE দশমিক
#অন্যথা
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
শূন্য NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  যদি (n <= 0) ফেরত;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  যদি (inc_in < 0)   in -= (n - 1) * inc_in;
  যদি (inc_out < 0) out -= (n - 1) * inc_out;

#যদি defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#নতুবা_যদি defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#অন্যথা
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#অন্যথা
শূন্য CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  যদি (n <= 0) ফেরত;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  যদি (inc_in < 0)   in -= (n - 1) * inc_in;
  যদি (inc_out < 0) out -= (n - 1) * inc_out;

#যদি defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#নতুবা_যদি defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#অন্যথা
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
