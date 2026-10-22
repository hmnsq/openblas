#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#যদি defined(DOUBLE_PREC)
#define FLOAT_TYPE দ্বিগুণ
#নোচেদ্যদি defined(SINGLE_PREC)
#define FLOAT_TYPE ভিন্ন
#অন্যথা
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
শূন্য NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   যদি (n <= 0) প্রত্যাগম;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   যদি (inc_in < 0)   in -= (n - 1) * inc_in;
   যদি (inc_out < 0) out -= (n - 1) * inc_out;

#যদি defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#নোচেদ্যদি defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#অন্যথা
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#অন্যথা
শূন্য CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  যদি (n <= 0) প্রত্যাগম;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  যদি (inc_in < 0)   in -= (n - 1) * inc_in;
  যদি (inc_out < 0) out -= (n - 1) * inc_out;

#যদি defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#নোচেদ্যদি defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
