#ضمن <stdio.h>
#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#إذا defined(DOUBLE_PREC)
#تعریف_کرو FLOAT_TYPE مزدوج
#وإلاإذا defined(SINGLE_PREC)
#تعریف_کرو FLOAT_TYPE عائم
#وإلا
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS
فراغ NAME(blasint *N, FLOAT_TYPE *وچ, blasint *INC_وچ, bfloat16 *باہرلا, blasint *INC_باہرلا){
   BLASLONG n    = *N;
   BLASLONG inc_وچ = *INC_وچ;
   BLASLONG inc_باہرلا = *INC_باہرلا;

   PRINT_DEBUG_NAME;

   إذا (n <= 0) أرجع;

   IDEBUG_شروع;
   FUNCTION_PROFILE_شروع();

   إذا (inc_وچ < 0)   وچ -= (n - 1) * inc_وچ;
   إذا (inc_باہرلا < 0) باہرلا -= (n - 1) * inc_باہرلا;

#إذا defined(DOUBLE_PREC)
   D_منزل_BF16_K(n, وچ, inc_وچ, باہرلا, inc_باہرلا);
#وإلاإذا defined(SINGLE_PREC)
   S_منزل_BF16_K(n, وچ, inc_وچ, باہرلا, inc_باہرلا);
#وإلا
#جے_انت

   FUNCTION_PROFILE_سرا(1, 2 * n, 2 * n);
   IDEBUG_سرا;
}
#وإلا
فراغ CNAME(blasint n, FLOAT_TYPE *وچ, blasint inc_وچ, bfloat16 *باہرلا, blasint inc_باہرلا){
  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع;

  IDEBUG_شروع;
  FUNCTION_PROFILE_شروع();

  إذا (inc_وچ < 0)   وچ -= (n - 1) * inc_وچ;
  إذا (inc_باہرلا < 0) باہرلا -= (n - 1) * inc_باہرلا;

#إذا defined(DOUBLE_PREC)
  D_منزل_BF16_K(n, وچ, inc_وچ, باہرلا, inc_باہرلا);
#وإلاإذا defined(SINGLE_PREC)
  S_منزل_BF16_K(n, وچ, inc_وچ, باہرلا, inc_باہرلا);
#جے_انت

  FUNCTION_PROFILE_سرا(1, 2 * n, 2 * n);
  IDEBUG_سرا;
}
#جے_انت
