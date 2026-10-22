#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ defined(DOUBLE_PREC)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ FLOAT_TYPE ਦੁੱਗਣਾ
#ਜਾਂਜੇ defined(SINGLE_PREC)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ FLOAT_TYPE ਦਸ਼ਮਲਵ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
ਖਾਲੀ NAME(blasint *N, FLOAT_TYPE *ਵਿੱਚ, blasint *INC_ਵਿੱਚ, bfloat16 *ਬਾਹਰਲਾ, blasint *INC_ਬਾਹਰਲਾ){
   BLASLONG n    = *N;
   BLASLONG inc_ਵਿੱਚ = *INC_ਵਿੱਚ;
   BLASLONG inc_ਬਾਹਰਲਾ = *INC_ਬਾਹਰਲਾ;

   PRINT_DEBUG_NAME;

   ਜੇ (n <= 0) ਵਾਪਸ;

   IDEBUG_ਸ਼ੁਰੂ;
   FUNCTION_PROFILE_ਸ਼ੁਰੂ();

   ਜੇ (inc_ਵਿੱਚ < 0)   ਵਿੱਚ -= (n - 1) * inc_ਵਿੱਚ;
   ਜੇ (inc_ਬਾਹਰਲਾ < 0) ਬਾਹਰਲਾ -= (n - 1) * inc_ਬਾਹਰਲਾ;

#ਜੇ defined(DOUBLE_PREC)
   D_ਮੰਜ਼ਿਲ_BF16_K(n, ਵਿੱਚ, inc_ਵਿੱਚ, ਬਾਹਰਲਾ, inc_ਬਾਹਰਲਾ);
#ਜਾਂਜੇ defined(SINGLE_PREC)
   S_ਮੰਜ਼ਿਲ_BF16_K(n, ਵਿੱਚ, inc_ਵਿੱਚ, ਬਾਹਰਲਾ, inc_ਬਾਹਰਲਾ);
#ਨਹੀਂਤਾਂ
#ਜੇ_ਅੰਤ

   FUNCTION_PROFILE_ਸਿਰਾ(1, 2 * n, 2 * n);
   IDEBUG_ਸਿਰਾ;
}
#ਨਹੀਂਤਾਂ
ਖਾਲੀ CNAME(blasint n, FLOAT_TYPE *ਵਿੱਚ, blasint inc_ਵਿੱਚ, bfloat16 *ਬਾਹਰਲਾ, blasint inc_ਬਾਹਰਲਾ){
  PRINT_DEBUG_CNAME;

  ਜੇ (n <= 0) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;
  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (inc_ਵਿੱਚ < 0)   ਵਿੱਚ -= (n - 1) * inc_ਵਿੱਚ;
  ਜੇ (inc_ਬਾਹਰਲਾ < 0) ਬਾਹਰਲਾ -= (n - 1) * inc_ਬਾਹਰਲਾ;

#ਜੇ defined(DOUBLE_PREC)
  D_ਮੰਜ਼ਿਲ_BF16_K(n, ਵਿੱਚ, inc_ਵਿੱਚ, ਬਾਹਰਲਾ, inc_ਬਾਹਰਲਾ);
#ਜਾਂਜੇ defined(SINGLE_PREC)
  S_ਮੰਜ਼ਿਲ_BF16_K(n, ਵਿੱਚ, inc_ਵਿੱਚ, ਬਾਹਰਲਾ, inc_ਬਾਹਰਲਾ);
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(1, 2 * n, 2 * n);
  IDEBUG_ਸਿਰਾ;
}
#ਜੇ_ਅੰਤ
