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
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QTRMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DTRMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "STRMV "
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  qtrmv_NUU, qtrmv_NUN, qtrmv_NLU, qtrmv_NLN,
  qtrmv_TUU, qtrmv_TUN, qtrmv_TLU, qtrmv_TLN,
#ਜਾਂਜੇ defined(DOUBLE)
  dtrmv_NUU, dtrmv_NUN, dtrmv_NLU, dtrmv_NLN,
  dtrmv_TUU, dtrmv_TUN, dtrmv_TLU, dtrmv_TLN,
#ਨਹੀਂਤਾਂ
  strmv_NUU, strmv_NUN, strmv_NLU, strmv_NLN,
  strmv_TUU, strmv_TUN, strmv_TLU, strmv_TLN,
#ਜੇ_ਅੰਤ
};

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  qtrmv_thread_NUU, qtrmv_thread_NUN, qtrmv_thread_NLU, qtrmv_thread_NLN,
  qtrmv_thread_TUU, qtrmv_thread_TUN, qtrmv_thread_TLU, qtrmv_thread_TLN,
#ਜਾਂਜੇ defined(DOUBLE)
  dtrmv_thread_NUU, dtrmv_thread_NUN, dtrmv_thread_NLU, dtrmv_thread_NLN,
  dtrmv_thread_TUU, dtrmv_thread_TUN, dtrmv_thread_TLU, dtrmv_thread_TLN,
#ਨਹੀਂਤਾਂ
  strmv_thread_NUU, strmv_thread_NUN, strmv_thread_NLU, strmv_thread_NLN,
  strmv_thread_TUU, strmv_thread_TUN, strmv_thread_TLU, strmv_thread_TLN,
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
  ਪੂਰਨ trans;
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
  ਜੇ (trans_arg == 'R') trans = 0;
  ਜੇ (trans_arg == 'C') trans = 1;

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
	   blasint n, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  ਪੂਰਨ trans, uplo, unit;
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
    ਜੇ (TransA == CblasConjNoTrans) trans = 0;
    ਜੇ (TransA == CblasConjTrans)   trans = 1;

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
    ਜੇ (TransA == CblasConjNoTrans) trans = 1;
    ਜੇ (TransA == CblasConjTrans)   trans = 0;

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

  ਜੇ (incx < 0 ) x -= (n - 1) * incx;

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  nthreads = ਅੰਕ_cpu_avail(2);
  
  ਜੇ (n < 50 ) nthreads = 1;
  ਜੇ (nthreads > 2 && n < 500) nthreads = 2;
	
  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

    (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(1, n * n / 2 + n,  n * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
