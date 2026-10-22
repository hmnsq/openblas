/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XTRMV "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZTRMV "
#وإلا
#تعریف_کن ERROR_NAME "CTRMV "
#پایان_اگر

ساكن صحيح (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#اگر_تعریف_شده XDOUBLE
  xtrmv_NUU, xtrmv_NUN, xtrmv_NLU, xtrmv_NLN,
  xtrmv_TUU, xtrmv_TUN, xtrmv_TLU, xtrmv_TLN,
  xtrmv_RUU, xtrmv_پیاپی, xtrmv_RLU, xtrmv_RLN,
  xtrmv_CUU, xtrmv_CUN, xtrmv_CLU, xtrmv_CLN,
#وإلاإذا defined(DOUBLE)
  ztrmv_NUU, ztrmv_NUN, ztrmv_NLU, ztrmv_NLN,
  ztrmv_TUU, ztrmv_TUN, ztrmv_TLU, ztrmv_TLN,
  ztrmv_RUU, ztrmv_پیاپی, ztrmv_RLU, ztrmv_RLN,
  ztrmv_CUU, ztrmv_CUN, ztrmv_CLU, ztrmv_CLN,
#وإلا
  ctrmv_NUU, ctrmv_NUN, ctrmv_NLU, ctrmv_NLN,
  ctrmv_TUU, ctrmv_TUN, ctrmv_TLU, ctrmv_TLN,
  ctrmv_RUU, ctrmv_پیاپی, ctrmv_RLU, ctrmv_RLN,
  ctrmv_CUU, ctrmv_CUN, ctrmv_CLU, ctrmv_CLN,
#پایان_اگر
};

#اگر_تعریف_شده SMP
ساكن صحيح (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#اگر_تعریف_شده XDOUBLE
  xtrmv_thread_NUU, xtrmv_thread_NUN, xtrmv_thread_NLU, xtrmv_thread_NLN,
  xtrmv_thread_TUU, xtrmv_thread_TUN, xtrmv_thread_TLU, xtrmv_thread_TLN,
  xtrmv_thread_RUU, xtrmv_thread_پیاپی, xtrmv_thread_RLU, xtrmv_thread_RLN,
  xtrmv_thread_CUU, xtrmv_thread_CUN, xtrmv_thread_CLU, xtrmv_thread_CLN,
#وإلاإذا defined(DOUBLE)
  ztrmv_thread_NUU, ztrmv_thread_NUN, ztrmv_thread_NLU, ztrmv_thread_NLN,
  ztrmv_thread_TUU, ztrmv_thread_TUN, ztrmv_thread_TLU, ztrmv_thread_TLN,
  ztrmv_thread_RUU, ztrmv_thread_پیاپی, ztrmv_thread_RLU, ztrmv_thread_RLN,
  ztrmv_thread_CUU, ztrmv_thread_CUN, ztrmv_thread_CLU, ztrmv_thread_CLN,
#وإلا
  ctrmv_thread_NUU, ctrmv_thread_NUN, ctrmv_thread_NLU, ctrmv_thread_NLN,
  ctrmv_thread_TUU, ctrmv_thread_TUN, ctrmv_thread_TLU, ctrmv_thread_TLN,
  ctrmv_thread_RUU, ctrmv_thread_پیاپی, ctrmv_thread_RLU, ctrmv_thread_RLN,
  ctrmv_thread_CUU, ctrmv_thread_CUN, ctrmv_thread_CLU, ctrmv_thread_CLN,
#پایان_اگر
};
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  صحيح uplo;
  صحيح unit;
  صحيح trans, میانگیر_اندازه;
  FLOAT *میانگیر;
#اگر_تعریف_شده SMP
  صحيح nthreads;
#پایان_اگر

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (diag_arg  == 'U') unit  = 0;
  إذا (diag_arg  == 'N') unit  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (incx == 0)          info =  8;
  إذا (lda  < MAX(1, n))   info =  6;
  إذا (n < 0)              info =  4;
  إذا (unit  < 0)          info =  3;
  إذا (trans < 0)          info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_DIAG Diag,
	   blasint n, فراغ  *va, blasint lda, فراغ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  صحيح trans, uplo, unit, میانگیر_اندازه;
  blasint info;
  FLOAT *میانگیر;
#اگر_تعریف_شده SMP
  صحيح nthreads;
#پایان_اگر

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Uplo == Cblas_فوقانی)         uplo  = 0;
    إذا (Uplo == Cblas_پایینی)         uplo  = 1;

    إذا (TransA == CblasNoTrans)     trans = 0;
    إذا (TransA == CblasTrans)       trans = 1;
    إذا (TransA == CblasConjNoTrans) trans = 2;
    إذا (TransA == CblasConjTrans)   trans = 3;

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    إذا (incx == 0)          info =  8;
    إذا (lda  < MAX(1, n))   info =  6;
    إذا (n < 0)              info =  4;
    إذا (unit  < 0)          info =  3;
    إذا (trans < 0)          info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتیب == Cblas_ردیف_Major) {
    إذا (Uplo == Cblas_فوقانی)         uplo  = 1;
    إذا (Uplo == Cblas_پایینی)         uplo  = 0;

    إذا (TransA == CblasNoTrans)     trans = 1;
    إذا (TransA == CblasTrans)       trans = 0;
    إذا (TransA == CblasConjNoTrans) trans = 3;
    إذا (TransA == CblasConjTrans)   trans = 2;

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    إذا (incx == 0)          info =  8;
    إذا (lda  < MAX(1, n))   info =  6;
    إذا (n < 0)              info =  4;
    إذا (unit  < 0)          info =  3;
    إذا (trans < 0)          info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#پایان_اگر

  إذا (n == 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;

#اگر_تعریف_شده SMP
  // Calibrated on a Xeon E5-2630
  إذا(1L * n * n > 36L * حجم(FLOAT) * حجم(FLOAT) * GEMM_MULTITHREAD_آستانه) {
    nthreads = رقم_cpu_avail(2);
    إذا(nthreads > 2 && 1L * n * n < 64L * حجم(FLOAT) * حجم(FLOAT) * GEMM_MULTITHREAD_آستانه)
      nthreads = 2;
  } وإلا
      nthreads = 1;

  إذا(nthreads > 1) {
    میانگیر_اندازه = n > 16 ? 0 : n * 4 + 40;
  }
  وإلا
#پایان_اگر
  {
    میانگیر_اندازه = ((n - 1) / DTB_ENTRIES) * 2 * DTB_ENTRIES + 32 / حجم(FLOAT);
    // It seems to be required for some K8 or Barcelona CPU
    میانگیر_اندازه += 8;
    إذا(incx != 1)
      میانگیر_اندازه += n * 2;
  }
  پشته_ALLOC(میانگیر_اندازه, FLOAT, میانگیر);

#اگر_تعریف_شده SMP
  إذا (nthreads == 1) {
#پایان_اگر

  (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, میانگیر);

#اگر_تعریف_شده SMP
  } وإلا {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, میانگیر, nthreads);

  }
#پایان_اگر

  پشته_FREE(میانگیر);

  FUNCTION_PROFILE_انتها(4, n * n / 2 + n,  n * n);

  IDEBUG_انتها;

  أرجع;
}
