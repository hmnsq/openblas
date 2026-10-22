#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि defined(DOUBLE_PREC)
#परिभाषय FLOAT_TYPE द्विगुण
#नोचेद्यदि defined(SINGLE_PREC)
#परिभाषय FLOAT_TYPE भिन्न
#अन्यथा
#यदि_अन्तः

OPENBLAS_EXPORT
#यदि_अपरिभाषितम् CBLAS
शून्य NAME(blasint *N, bfloat16 *अन्तर्, blasint *INC_अन्तर्, FLOAT_TYPE *बहिर्गमनम्, blasint *INC_बहिर्गमनम्){
  BLASLONG n    = *N;
  BLASLONG inc_अन्तर् = *INC_अन्तर्;
  BLASLONG inc_बहिर्गमनम् = *INC_बहिर्गमनम्;

  PRINT_DEBUG_NAME;

  यदि (n <= 0) प्रत्यागम;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  यदि (inc_अन्तर् < 0)   अन्तर् -= (n - 1) * inc_अन्तर्;
  यदि (inc_बहिर्गमनम् < 0) बहिर्गमनम् -= (n - 1) * inc_बहिर्गमनम्;

#यदि defined(DOUBLE_PREC)
  D_BF16_गन्तव्यम्_K(n, अन्तर्, inc_अन्तर्, बहिर्गमनम्, inc_बहिर्गमनम्);
#नोचेद्यदि defined(SINGLE_PREC)
  S_BF16_गन्तव्यम्_K(n, अन्तर्, inc_अन्तर्, बहिर्गमनम्, inc_बहिर्गमनम्);
#अन्यथा
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(1, 2 * n, 2 * n);
  IDEBUG_अन्तः;
}
#अन्यथा
शून्य CNAME(blasint n, bfloat16 * अन्तर्, blasint inc_अन्तर्, FLOAT_TYPE * बहिर्गमनम्, blasint inc_बहिर्गमनम्){
  PRINT_DEBUG_CNAME;

  यदि (n <= 0) प्रत्यागम;

  IDEBUG_START;
  FUNCTION_PROFILE_START();

  यदि (inc_अन्तर् < 0)   अन्तर् -= (n - 1) * inc_अन्तर्;
  यदि (inc_बहिर्गमनम् < 0) बहिर्गमनम् -= (n - 1) * inc_बहिर्गमनम्;

#यदि defined(DOUBLE_PREC)
  D_BF16_गन्तव्यम्_K(n, अन्तर्, inc_अन्तर्, बहिर्गमनम्, inc_बहिर्गमनम्);
#नोचेद्यदि defined(SINGLE_PREC)
  S_BF16_गन्तव्यम्_K(n, अन्तर्, inc_अन्तर्, बहिर्गमनम्, inc_बहिर्गमनम्);
#अन्यथा
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(1, 2 * n, 2 * n);
  IDEBUG_अन्तः;
}
#यदि_अन्तः
