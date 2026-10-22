#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
ದಶಮಾಂಶ NAME(blasint *N, bfloat16 *x, blasint *INCX, bfloat16 *y, blasint *INCY){
   BLASLONG n    = *N;
   BLASLONG incx = *INCX;
   BLASLONG incy = *INCY;
   ದಶಮಾಂಶ ret;
   PRINT_DEBUG_NAME;

   ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ 0.;

   IDEBUG_ಆರಂಭ;
   FUNCTION_PROFILE_ಆರಂಭ();

   ಆದರೆ (incx < 0) x -= (n - 1) * incx;
   ಆದರೆ (incy < 0) y -= (n - 1) * incy;
   ret = BF16_DOT_K(n, x, incx, y, incy);

   FUNCTION_PROFILE_ಕೊನೆ(1, 2 * n, 2 * n);
   IDEBUG_ಕೊನೆ;

   ಹಿಂತಿರುಗಿ ret;
 }

#ಇಲ್ಲದಿದ್ದರೆ

ದಶಮಾಂಶ CNAME(blasint n, bfloat16 *x, blasint incx, bfloat16 *y, blasint incy){

  ದಶಮಾಂಶ ret;
  PRINT_DEBUG_CNAME;

  ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ 0.;

  IDEBUG_ಆರಂಭ;
  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0) x -= (n - 1) * incx;
  ಆದರೆ (incy < 0) y -= (n - 1) * incy;
  ret = BF16_DOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_ಕೊನೆ(1, 2 * n, 2 * n);
  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ ret;
}

#ಆದರೆ_ಅಂತ್ಯ
