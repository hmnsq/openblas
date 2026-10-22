#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر defined(DOUBLE_PREC)
#تعریف_کرو FLOAT_TYPE دوہرا
#وگرنہ_اگر defined(SINGLE_PREC)
#تعریف_کرو FLOAT_TYPE اعشاری
#ورنہ
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS
خالی NAME(blasint *N, FLOAT_TYPE *اندر, blasint *INC_اندر, bfloat16 *برآمد, blasint *INC_برآمد){
   BLASLONG n    = *N;
   BLASLONG inc_اندر = *INC_اندر;
   BLASLONG inc_برآمد = *INC_برآمد;

   PRINT_DEBUG_NAME;

   اگر (n <= 0) واپس;

   IDEBUG_آغاز;
   FUNCTION_PROFILE_آغاز();

   اگر (inc_اندر < 0)   اندر -= (n - 1) * inc_اندر;
   اگر (inc_برآمد < 0) برآمد -= (n - 1) * inc_برآمد;

#اگر defined(DOUBLE_PREC)
   D_مقصد_BF16_K(n, اندر, inc_اندر, برآمد, inc_برآمد);
#وگرنہ_اگر defined(SINGLE_PREC)
   S_مقصد_BF16_K(n, اندر, inc_اندر, برآمد, inc_برآمد);
#ورنہ
#اگر_ختم

   FUNCTION_PROFILE_انتہا(1, 2 * n, 2 * n);
   IDEBUG_انتہا;
}
#ورنہ
خالی CNAME(blasint n, FLOAT_TYPE *اندر, blasint inc_اندر, bfloat16 *برآمد, blasint inc_برآمد){
  PRINT_DEBUG_CNAME;

  اگر (n <= 0) واپس;

  IDEBUG_آغاز;
  FUNCTION_PROFILE_آغاز();

  اگر (inc_اندر < 0)   اندر -= (n - 1) * inc_اندر;
  اگر (inc_برآمد < 0) برآمد -= (n - 1) * inc_برآمد;

#اگر defined(DOUBLE_PREC)
  D_مقصد_BF16_K(n, اندر, inc_اندر, برآمد, inc_برآمد);
#وگرنہ_اگر defined(SINGLE_PREC)
  S_مقصد_BF16_K(n, اندر, inc_اندر, برآمد, inc_برآمد);
#اگر_ختم

  FUNCTION_PROFILE_انتہا(1, 2 * n, 2 * n);
  IDEBUG_انتہا;
}
#اگر_ختم
