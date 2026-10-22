#ضمن <stdio.h>
#ضمن "common.h"
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#إذا defined(DOUBLE_PREC)
#وصف_ڪريو FLOAT_TYPE مزدوج
#وإلاإذا defined(SINGLE_PREC)
#وصف_ڪريو FLOAT_TYPE عائم
#وإلا
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS
فراغ NAME(blasint *N, FLOAT_TYPE *اندريون, blasint *INC_اندريون, bfloat16 *ٻاهريون, blasint *INC_ٻاهريون){
   BLASLONG n    = *N;
   BLASLONG inc_اندريون = *INC_اندريون;
   BLASLONG inc_ٻاهريون = *INC_ٻاهريون;

   PRINT_DEBUG_NAME;

   إذا (n <= 0) أرجع;

   IDEBUG_شروعات;
   FUNCTION_PROFILE_شروعات();

   إذا (inc_اندريون < 0)   اندريون -= (n - 1) * inc_اندريون;
   إذا (inc_ٻاهريون < 0) ٻاهريون -= (n - 1) * inc_ٻاهريون;

#إذا defined(DOUBLE_PREC)
   D_منزل_BF16_K(n, اندريون, inc_اندريون, ٻاهريون, inc_ٻاهريون);
#وإلاإذا defined(SINGLE_PREC)
   S_منزل_BF16_K(n, اندريون, inc_اندريون, ٻاهريون, inc_ٻاهريون);
#وإلا
#جيڪڏهن_پورو

   FUNCTION_PROFILE_پڇاڙي(1, 2 * n, 2 * n);
   IDEBUG_پڇاڙي;
}
#وإلا
فراغ CNAME(blasint n, FLOAT_TYPE *اندريون, blasint inc_اندريون, bfloat16 *ٻاهريون, blasint inc_ٻاهريون){
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع;

  IDEBUG_شروعات;
  FUNCTION_PROFILE_شروعات();

  إذا (inc_اندريون < 0)   اندريون -= (n - 1) * inc_اندريون;
  إذا (inc_ٻاهريون < 0) ٻاهريون -= (n - 1) * inc_ٻاهريون;

#إذا defined(DOUBLE_PREC)
  D_منزل_BF16_K(n, اندريون, inc_اندريون, ٻاهريون, inc_ٻاهريون);
#وإلاإذا defined(SINGLE_PREC)
  S_منزل_BF16_K(n, اندريون, inc_اندريون, ٻاهريون, inc_ٻاهريون);
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(1, 2 * n, 2 * n);
  IDEBUG_پڇاڙي;
}
#جيڪڏهن_پورو
