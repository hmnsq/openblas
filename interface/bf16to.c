#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#जर defined(DOUBLE_PREC)
#define FLOAT_TYPE द्विगुण
#किंवाजर defined(SINGLE_PREC)
#define FLOAT_TYPE दशांश
#अन्यथा
#endif

OPENBLAS_EXPORT
#ifndef CBLAS
शून्य NAME(blasint *N, bfloat16 *in, blasint *INC_IN, FLOAT_TYPE *out, blasint *INC_OUT){
  BLASLONG n    = *N;
  BLASLONG inc_in = *INC_IN;
  BLASLONG inc_out = *INC_OUT;

  PRINT_DEBUG_NAME;

  जर (n <= 0) परत;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  जर (inc_in < 0)   in -= (n - 1) * inc_in;
  जर (inc_out < 0) out -= (n - 1) * inc_out;

#जर defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#किंवाजर defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#अन्यथा
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#अन्यथा
शून्य CNAME(blasint n, bfloat16 * in, blasint inc_in, FLOAT_TYPE * out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  जर (n <= 0) परत;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  जर (inc_in < 0)   in -= (n - 1) * inc_in;
  जर (inc_out < 0) out -= (n - 1) * inc_out;

#जर defined(DOUBLE_PREC)
  D_BF16_TO_K(n, in, inc_in, out, inc_out);
#किंवाजर defined(SINGLE_PREC)
  S_BF16_TO_K(n, in, inc_in, out, inc_out);
#अन्यथा
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#endif
