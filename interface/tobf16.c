#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर defined(DOUBLE_PREC)
#घोषणा FLOAT_TYPE द्विगुन
#या_अगर defined(SINGLE_PREC)
#घोषणा FLOAT_TYPE भग्न
#अन्यथा
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS
व्योम NAME(blasint *N, FLOAT_TYPE *in, blasint *INC_IN, bfloat16 *out, blasint *INC_OUT){
   BLASLONG n    = *N;
   BLASLONG inc_in = *INC_IN;
   BLASLONG inc_out = *INC_OUT;

   PRINT_DEBUG_NAME;

   अगर (n <= 0) वापस;

   IDEBUG_START;
   FUNCTION_PROFILE_START();

   अगर (inc_in < 0)   in -= (n - 1) * inc_in;
   अगर (inc_out < 0) out -= (n - 1) * inc_out;

#अगर defined(DOUBLE_PREC)
   D_TO_BF16_K(n, in, inc_in, out, inc_out);
#या_अगर defined(SINGLE_PREC)
   S_TO_BF16_K(n, in, inc_in, out, inc_out);
#अन्यथा
#पूर्ण_अगर

   FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
   IDEBUG_END;
}
#अन्यथा
व्योम CNAME(blasint n, FLOAT_TYPE *in, blasint inc_in, bfloat16 *out, blasint inc_out){
  PRINT_DEBUG_CNAME;

  अगर (n <= 0) वापस;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  अगर (inc_in < 0)   in -= (n - 1) * inc_in;
  अगर (inc_out < 0) out -= (n - 1) * inc_out;

#अगर defined(DOUBLE_PREC)
  D_TO_BF16_K(n, in, inc_in, out, inc_out);
#या_अगर defined(SINGLE_PREC)
  S_TO_BF16_K(n, in, inc_in, out, inc_out);
#पूर्ण_अगर

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);
  IDEBUG_END;
}
#पूर्ण_अगर
