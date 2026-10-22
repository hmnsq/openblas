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
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XTPSV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZTPSV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CTPSV "
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, ਖਾਲੀ *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xtpsv_NUU, xtpsv_NUN, xtpsv_NLU, xtpsv_NLN,
  xtpsv_TUU, xtpsv_TUN, xtpsv_TLU, xtpsv_TLN,
  xtpsv_RUU, xtpsv_ਲੜੀਵਾਰ, xtpsv_RLU, xtpsv_RLN,
  xtpsv_CUU, xtpsv_CUN, xtpsv_CLU, xtpsv_CLN,
#ਜਾਂਜੇ defined(DOUBLE)
  ztpsv_NUU, ztpsv_NUN, ztpsv_NLU, ztpsv_NLN,
  ztpsv_TUU, ztpsv_TUN, ztpsv_TLU, ztpsv_TLN,
  ztpsv_RUU, ztpsv_ਲੜੀਵਾਰ, ztpsv_RLU, ztpsv_RLN,
  ztpsv_CUU, ztpsv_CUN, ztpsv_CLU, ztpsv_CLN,
#ਨਹੀਂਤਾਂ
  ctpsv_NUU, ctpsv_NUN, ctpsv_NLU, ctpsv_NLN,
  ctpsv_TUU, ctpsv_TUN, ctpsv_TLU, ctpsv_TLN,
  ctpsv_RUU, ctpsv_ਲੜੀਵਾਰ, ctpsv_RLU, ctpsv_RLN,
  ctpsv_CUU, ctpsv_CUN, ctpsv_CLU, ctpsv_CLN,
#ਜੇ_ਅੰਤ
};

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, ਅੱਖਰ *TRANS, ਅੱਖਰ *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  ਅੱਖਰ uplo_arg  = *UPLO;
  ਅੱਖਰ trans_arg = *TRANS;
  ਅੱਖਰ diag_arg  = *DIAG;

  blasint n    = *N;
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

  ਜੇ (incx == 0)          info =  7;
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
	   blasint n, ਖਾਲੀ  *va, ਖਾਲੀ  *vx, blasint incx) {

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

    ਜੇ (incx == 0)          info =  7;
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

    ਜੇ (incx == 0)          info =  7;
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

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, ਬਫਰ);

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(4, n * n / 2 + n, n * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
