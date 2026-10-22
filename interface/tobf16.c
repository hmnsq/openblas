#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ଯଦି defined(DOUBLE_PREC)
#define FLOAT_TYPE ଦ୍ୱିଗୁଣ
#ନଚେତଯଦି defined(SINGLE_PREC)
#define FLOAT_TYPE ଦଶମିକ
#ଅନ୍ୟଥା
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
ଶୂନ୍ୟ NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   ଯଦି (n <= 0) ଫେରନ୍ତୁ;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   ଯଦି (inc_in < 0)   in -= (n - 1) * inc_in;
   ଯଦି (inc_out < 0) out -= (n - 1) * inc_out;

#ଯଦି defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#ନଚେତଯଦି defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#ଅନ୍ୟଥା
#endif

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#ଅନ୍ୟଥା
ଶୂନ୍ୟ CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  ଯଦି (n <= 0) ଫେରନ୍ତୁ;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  ଯଦି (inc_in < 0)   in -= (n - 1) * inc_in;
  ଯଦି (inc_out < 0) out -= (n - 1) * inc_out;

#ଯଦି defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#ନଚେତଯଦି defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
