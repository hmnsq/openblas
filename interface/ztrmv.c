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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <ctype.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XTRMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZTRMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CTRMV "
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xtrmv_NUU, xtrmv_NUN, xtrmv_NLU, xtrmv_NLN,
  xtrmv_TUU, xtrmv_TUN, xtrmv_TLU, xtrmv_TLN,
  xtrmv_RUU, xtrmv_ਲੜੀਵਾਰ, xtrmv_RLU, xtrmv_RLN,
  xtrmv_CUU, xtrmv_CUN, xtrmv_CLU, xtrmv_CLN,
#ਜਾਂਜੇ defined(DOUBLE)
  ztrmv_NUU, ztrmv_NUN, ztrmv_NLU, ztrmv_NLN,
  ztrmv_TUU, ztrmv_TUN, ztrmv_TLU, ztrmv_TLN,
  ztrmv_RUU, ztrmv_ਲੜੀਵਾਰ, ztrmv_RLU, ztrmv_RLN,
  ztrmv_CUU, ztrmv_CUN, ztrmv_CLU, ztrmv_CLN,
#ਨਹੀਂਤਾਂ
  ctrmv_NUU, ctrmv_NUN, ctrmv_NLU, ctrmv_NLN,
  ctrmv_TUU, ctrmv_TUN, ctrmv_TLU, ctrmv_TLN,
  ctrmv_RUU, ctrmv_ਲੜੀਵਾਰ, ctrmv_RLU, ctrmv_RLN,
  ctrmv_CUU, ctrmv_CUN, ctrmv_CLU, ctrmv_CLN,
#ਜੇ_ਅੰਤ
};

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xtrmv_thread_NUU, xtrmv_thread_NUN, xtrmv_thread_NLU, xtrmv_thread_NLN,
  xtrmv_thread_TUU, xtrmv_thread_TUN, xtrmv_thread_TLU, xtrmv_thread_TLN,
  xtrmv_thread_RUU, xtrmv_thread_ਲੜੀਵਾਰ, xtrmv_thread_RLU, xtrmv_thread_RLN,
  xtrmv_thread_CUU, xtrmv_thread_CUN, xtrmv_thread_CLU, xtrmv_thread_CLN,
#ਜਾਂਜੇ defined(DOUBLE)
  ztrmv_thread_NUU, ztrmv_thread_NUN, ztrmv_thread_NLU, ztrmv_thread_NLN,
  ztrmv_thread_TUU, ztrmv_thread_TUN, ztrmv_thread_TLU, ztrmv_thread_TLN,
  ztrmv_thread_RUU, ztrmv_thread_ਲੜੀਵਾਰ, ztrmv_thread_RLU, ztrmv_thread_RLN,
  ztrmv_thread_CUU, ztrmv_thread_CUN, ztrmv_thread_CLU, ztrmv_thread_CLN,
#ਨਹੀਂਤਾਂ
  ctrmv_thread_NUU, ctrmv_thread_NUN, ctrmv_thread_NLU, ctrmv_thread_NLN,
  ctrmv_thread_TUU, ctrmv_thread_TUN, ctrmv_thread_TLU, ctrmv_thread_TLN,
  ctrmv_thread_RUU, ctrmv_thread_ਲੜੀਵਾਰ, ctrmv_thread_RLU, ctrmv_thread_RLN,
  ctrmv_thread_CUU, ctrmv_thread_CUN, ctrmv_thread_CLU, ctrmv_thread_CLN,
#ਜੇ_ਅੰਤ
};
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, ਅੱਖਰ *TRANS, ਅੱਖਰ *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  ਅੱਖਰ uplo_arg  = *UPLO;
  ਅੱਖਰ trans_arg = *TRANS;
  ਅੱਖਰ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  ਪੂਰਨ uplo;
  ਪੂਰਨ unit;
  ਪੂਰਨ trans, ਬਫਰ_ਮਾਪ;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ਜੇ (trans_arg == 'N') trans = 0;
  ਜੇ (trans_arg == 'T') trans = 1;
  ਜੇ (trans_arg == 'R') trans = 2;
  ਜੇ (trans_arg == 'C') trans = 3;

  ਜੇ (diag_arg  == 'U') unit  = 0;
  ਜੇ (diag_arg  == 'N') unit  = 1;

  ਜੇ (uplo_arg  == 'U') uplo  = 0;
  ਜੇ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ਜੇ (incx == 0)          info =  8;
  ਜੇ (lda  < MAX(1, n))   info =  6;
  ਜੇ (n < 0)              info =  4;
  ਜੇ (unit  < 0)          info =  3;
  ਜੇ (trans < 0)          info =  2;
  ਜੇ (uplo  < 0)          info =  1;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_UPLO Uplo,
	   ਗਿਣਤੀ CBLAS_TRANSPOSE TransA, ਗਿਣਤੀ CBLAS_DIAG Diag,
	   blasint n, ਖਾਲੀ  *va, blasint lda, ਖਾਲੀ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  ਪੂਰਨ trans, uplo, unit, ਬਫਰ_ਮਾਪ;
  blasint info;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ)         uplo  = 0;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ)         uplo  = 1;

    ਜੇ (TransA == CblasNoTrans)     trans = 0;
    ਜੇ (TransA == CblasTrans)       trans = 1;
    ਜੇ (TransA == CblasConjNoTrans) trans = 2;
    ਜੇ (TransA == CblasConjTrans)   trans = 3;

    ਜੇ (Diag == CblasUnit)          unit  = 0;
    ਜੇ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ਜੇ (incx == 0)          info =  8;
    ਜੇ (lda  < MAX(1, n))   info =  6;
    ਜੇ (n < 0)              info =  4;
    ਜੇ (unit  < 0)          info =  3;
    ਜੇ (trans < 0)          info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ)         uplo  = 1;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ)         uplo  = 0;

    ਜੇ (TransA == CblasNoTrans)     trans = 1;
    ਜੇ (TransA == CblasTrans)       trans = 0;
    ਜੇ (TransA == CblasConjNoTrans) trans = 3;
    ਜੇ (TransA == CblasConjTrans)   trans = 2;

    ਜੇ (Diag == CblasUnit)          unit  = 0;
    ਜੇ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ਜੇ (incx == 0)          info =  8;
    ਜੇ (lda  < MAX(1, n))   info =  6;
    ਜੇ (n < 0)              info =  4;
    ਜੇ (unit  < 0)          info =  3;
    ਜੇ (trans < 0)          info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  ਜੇ (n == 0) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0 ) x -= (n - 1) * incx * 2;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  // Calibrated on a Xeon E5-2630
  ਜੇ(1L * n * n > 36L * ਆਕਾਰ(FLOAT) * ਆਕਾਰ(FLOAT) * GEMM_MULTITHREAD_ਹੱਦ) {
    nthreads = ਅੰਕ_cpu_avail(2);
    ਜੇ(nthreads > 2 && 1L * n * n < 64L * ਆਕਾਰ(FLOAT) * ਆਕਾਰ(FLOAT) * GEMM_MULTITHREAD_ਹੱਦ)
      nthreads = 2;
  } ਨਹੀਂਤਾਂ
      nthreads = 1;

  ਜੇ(nthreads > 1) {
    ਬਫਰ_ਮਾਪ = n > 16 ? 0 : n * 4 + 40;
  }
  ਨਹੀਂਤਾਂ
#ਜੇ_ਅੰਤ
  {
    ਬਫਰ_ਮਾਪ = ((n - 1) / DTB_ENTRIES) * 2 * DTB_ENTRIES + 32 / ਆਕਾਰ(FLOAT);
    // It seems to be required for some K8 or Barcelona CPU
    ਬਫਰ_ਮਾਪ += 8;
    ਜੇ(incx != 1)
      ਬਫਰ_ਮਾਪ += n * 2;
  }
  ਢੇਰ_ALLOC(ਬਫਰ_ਮਾਪ, FLOAT, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  ਢੇਰ_FREE(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(4, n * n / 2 + n,  n * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
