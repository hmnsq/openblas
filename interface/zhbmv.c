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
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XHBMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZHBMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CHBMV "
#ਜੇ_ਅੰਤ

ਸਥਿਰ  ਪੂਰਨ (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ਖਾਲੀ *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#ਜਾਂਜੇ defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#ਨਹੀਂਤਾਂ
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#ਜੇ_ਅੰਤ
};

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPBUG
ਸਥਿਰ  ਪੂਰਨ (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#ਜਾਂਜੇ defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#ਨਹੀਂਤਾਂ
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#ਜੇ_ਅੰਤ
};
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ਅੱਖਰ uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  ਪੂਰਨ uplo;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPBUG
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ਜੇ (uplo_arg  == 'U') uplo  = 0;
  ਜੇ (uplo_arg  == 'L') uplo  = 1;
  ਜੇ (uplo_arg  == 'V') uplo  = 2;
  ਜੇ (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  ਜੇ (incy == 0)          info = 11;
  ਜੇ (incx == 0)          info =  8;
  ਜੇ (lda  < k + 1)       info =  6;
  ਜੇ (k < 0)              info =  3;
  ਜੇ (n < 0)              info =  2;
  ਜੇ (uplo  < 0)          info =  1;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ,
	   ਗਿਣਤੀ CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   ਖਾਲੀ *VALPHA,
	   ਖਾਲੀ  *va, blasint lda,
	   ਖਾਲੀ  *vx, blasint incx,
	   ਖਾਲੀ *VBETA,
	   ਖਾਲੀ  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  FLOAT *ਬਫਰ;
  ਪੂਰਨ uplo;
  blasint info;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPBUG
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ)         uplo  = 0;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ)         uplo  = 1;

    info = -1;

    ਜੇ (incy == 0)          info = 11;
    ਜੇ (incx == 0)          info =  8;
    ਜੇ (lda  < k + 1)       info =  6;
    ਜੇ (k < 0)              info =  3;
    ਜੇ (n < 0)              info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ)         uplo  = 3;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ)         uplo  = 2;

    info = -1;

    ਜੇ (incy == 0)          info = 11;
    ਜੇ (incx == 0)          info =  8;
    ਜੇ (lda  < k + 1)       info =  6;
    ਜੇ (k < 0)              info =  3;
    ਜੇ (n < 0)              info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  ਜੇ (n == 0) ਵਾਪਸ;

  ਜੇ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ਨਿਰਾ, 0, ਨਿਰਾ, 0);

  ਜੇ ((alpha_r == ZERO) && (alpha_i == ZERO)) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  ਜੇ (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPBUG
  nthreads = ਅੰਕ_cpu_avail(2);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPBUG
  } ਨਹੀਂਤਾਂ {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(4, n * k / 2 + n,  n * k);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
