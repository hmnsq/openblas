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
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XTBSV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZTBSV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CTBSV "
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ਖਾਲੀ *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xtbsv_NUU, xtbsv_NUN, xtbsv_NLU, xtbsv_NLN,
  xtbsv_TUU, xtbsv_TUN, xtbsv_TLU, xtbsv_TLN,
  xtbsv_RUU, xtbsv_ਲੜੀਵਾਰ, xtbsv_RLU, xtbsv_RLN,
  xtbsv_CUU, xtbsv_CUN, xtbsv_CLU, xtbsv_CLN,
#ਜਾਂਜੇ defined(DOUBLE)
  ztbsv_NUU, ztbsv_NUN, ztbsv_NLU, ztbsv_NLN,
  ztbsv_TUU, ztbsv_TUN, ztbsv_TLU, ztbsv_TLN,
  ztbsv_RUU, ztbsv_ਲੜੀਵਾਰ, ztbsv_RLU, ztbsv_RLN,
  ztbsv_CUU, ztbsv_CUN, ztbsv_CLU, ztbsv_CLN,
#ਨਹੀਂਤਾਂ
  ctbsv_NUU, ctbsv_NUN, ctbsv_NLU, ctbsv_NLN,
  ctbsv_TUU, ctbsv_TUN, ctbsv_TLU, ctbsv_TLN,
  ctbsv_RUU, ctbsv_ਲੜੀਵਾਰ, ctbsv_RLU, ctbsv_RLN,
  ctbsv_CUU, ctbsv_CUN, ctbsv_CLU, ctbsv_CLN,
#ਜੇ_ਅੰਤ
};

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, ਅੱਖਰ *TRANS, ਅੱਖਰ *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  ਅੱਖਰ uplo_arg  = *UPLO;
  ਅੱਖਰ trans_arg = *TRANS;
  ਅੱਖਰ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  ਪੂਰਨ uplo;
  ਪੂਰਨ unit;
  ਪੂਰਨ trans;
  FLOAT *ਬਫਰ;

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

  ਜੇ (incx == 0)          info =  9;
  ਜੇ (lda < k + 1)        info =  7;
  ਜੇ (k < 0)              info =  5;
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
	   blasint n, blasint k, ਖਾਲੀ  *va, blasint lda, ਖਾਲੀ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  ਪੂਰਨ trans, uplo, unit;
  blasint info;
  FLOAT *ਬਫਰ;

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

    ਜੇ (incx == 0)          info =  9;
    ਜੇ (lda < k + 1)        info =  7;
    ਜੇ (k < 0)              info =  5;
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

    ਜੇ (incx == 0)          info =  9;
    ਜੇ (lda < k + 1)        info =  7;
    ਜੇ (k < 0)              info =  5;
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

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, ਬਫਰ);

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(4, n * k / 2 + n, n * k);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
