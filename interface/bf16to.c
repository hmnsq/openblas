#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(DOUBLE_PREC)
#ವ್ಯಾಖ್ಯಾನಿಸು FLOAT_TYPE ದ್ವಿಗುಣ
#ಅಥವಾಆದರೆ defined(SINGLE_PREC)
#ವ್ಯಾಖ್ಯಾನಿಸು FLOAT_TYPE ದಶಮಾಂಶ
#ಇಲ್ಲದಿದ್ದರೆ
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
ಶೂನ್ಯ NAME(blasint *N, bfloat16 *ಒಳ, blasint *INC_ಒಳ, FLOAT_TYPE *ಹೊರಹೋಗು, blasint *INC_ಹೊರಹೋಗು){
  BLASLONG n    = *N;
  BLASLONG inc_ಒಳ = *INC_ಒಳ;
  BLASLONG inc_ಹೊರಹೋಗು = *INC_ಹೊರಹೋಗು;

  PRINT_DEBUG_NAME;

  ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;
  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (inc_ಒಳ < 0)   ಒಳ -= (n - 1) * inc_ಒಳ;
  ಆದರೆ (inc_ಹೊರಹೋಗು < 0) ಹೊರಹೋಗು -= (n - 1) * inc_ಹೊರಹೋಗು;

#ಆದರೆ defined(DOUBLE_PREC)
  D_BF16_ಗಮ್ಯ_K(n, ಒಳ, inc_ಒಳ, ಹೊರಹೋಗು, inc_ಹೊರಹೋಗು);
#ಅಥವಾಆದರೆ defined(SINGLE_PREC)
  S_BF16_ಗಮ್ಯ_K(n, ಒಳ, inc_ಒಳ, ಹೊರಹೋಗು, inc_ಹೊರಹೋಗು);
#ಇಲ್ಲದಿದ್ದರೆ
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(1, 2 * n, 2 * n);
  IDEBUG_ಕೊನೆ;
}
#ಇಲ್ಲದಿದ್ದರೆ
ಶೂನ್ಯ CNAME(blasint n, bfloat16 * ಒಳ, blasint inc_ಒಳ, FLOAT_TYPE * ಹೊರಹೋಗು, blasint inc_ಹೊರಹೋಗು){
  PRINT_DEBUG_CNAME;

  ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;
  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (inc_ಒಳ < 0)   ಒಳ -= (n - 1) * inc_ಒಳ;
  ಆದರೆ (inc_ಹೊರಹೋಗು < 0) ಹೊರಹೋಗು -= (n - 1) * inc_ಹೊರಹೋಗು;

#ಆದರೆ defined(DOUBLE_PREC)
  D_BF16_ಗಮ್ಯ_K(n, ಒಳ, inc_ಒಳ, ಹೊರಹೋಗು, inc_ಹೊರಹೋಗು);
#ಅಥವಾಆದರೆ defined(SINGLE_PREC)
  S_BF16_ಗಮ್ಯ_K(n, ಒಳ, inc_ಒಳ, ಹೊರಹೋಗು, inc_ಹೊರಹೋಗು);
#ಇಲ್ಲದಿದ್ದರೆ
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(1, 2 * n, 2 * n);
  IDEBUG_ಕೊನೆ;
}
#ಆದರೆ_ಅಂತ್ಯ
