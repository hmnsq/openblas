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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <ctype.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHBMV "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZHBMV "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CHBMV "
#endif

സ്ഥിതം  പൂർണ്ണസംഖ്യ (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ശൂന്യം *) = {
#ifdef XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#അല്ലെങ്കിൽ
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#endif
};

#ifdef SMPBUG
സ്ഥിതം  പൂർണ്ണസംഖ്യ (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#അല്ലെങ്കിൽ
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  അക്ഷരം uplo_arg = *UPLO;
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
  പൂർണ്ണസംഖ്യ uplo;
  FLOAT *buffer;
#ifdef SMPBUG
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  എങ്കിൽ (uplo_arg  == 'U') uplo  = 0;
  എങ്കിൽ (uplo_arg  == 'L') uplo  = 1;
  എങ്കിൽ (uplo_arg  == 'V') uplo  = 2;
  എങ്കിൽ (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  എങ്കിൽ (incy == 0)          info = 11;
  എങ്കിൽ (incx == 0)          info =  8;
  എങ്കിൽ (lda  < k + 1)       info =  6;
  എങ്കിൽ (k < 0)              info =  3;
  എങ്കിൽ (n < 0)              info =  2;
  എങ്കിൽ (uplo  < 0)          info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order,
	   എണ്ണൽ CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   ശൂന്യം *VALPHA,
	   ശൂന്യം  *va, blasint lda,
	   ശൂന്യം  *vx, blasint incx,
	   ശൂന്യം *VBETA,
	   ശൂന്യം  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  FLOAT *buffer;
  പൂർണ്ണസംഖ്യ uplo;
  blasint info;
#ifdef SMPBUG
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (Uplo == CblasUpper)         uplo  = 0;
    എങ്കിൽ (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    എങ്കിൽ (incy == 0)          info = 11;
    എങ്കിൽ (incx == 0)          info =  8;
    എങ്കിൽ (lda  < k + 1)       info =  6;
    എങ്കിൽ (k < 0)              info =  3;
    എങ്കിൽ (n < 0)              info =  2;
    എങ്കിൽ (uplo  < 0)          info =  1;
  }

  എങ്കിൽ (order == CblasRowMajor) {
    എങ്കിൽ (Uplo == CblasUpper)         uplo  = 3;
    എങ്കിൽ (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    എങ്കിൽ (incy == 0)          info = 11;
    എങ്കിൽ (incx == 0)          info =  8;
    എങ്കിൽ (lda  < k + 1)       info =  6;
    എങ്കിൽ (k < 0)              info =  3;
    എങ്കിൽ (n < 0)              info =  2;
    എങ്കിൽ (uplo  < 0)          info =  1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#endif

  എങ്കിൽ (n == 0) തിരികെ;

  എങ്കിൽ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ശൂന്യ, 0, ശൂന്യ, 0);

  എങ്കിൽ ((alpha_r == ZERO) && (alpha_i == ZERO)) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  എങ്കിൽ (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  എങ്കിൽ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#ifdef SMPBUG
  nthreads = അക്കം_cpu_avail(2);

  എങ്കിൽ (nthreads == 1) {
#endif

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMPBUG
  } അല്ലെങ്കിൽ {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n,  n * k);

  IDEBUG_END;

  തിരികെ;
}
