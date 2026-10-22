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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <ctype.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHBMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHBMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHBMV "
#endif

ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#ಅನ್ಯಥಾ
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#endif
};

#ifdef SMPBUG
ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#ಅನ್ಯಥಾ
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ವರ್ಣ uplo_arg = *UPLO;
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
  ಪೂರ್ಣಾಙ್ಕ uplo;
  FLOAT *buffer;
#ifdef SMPBUG
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;
  ಯದಿ (uplo_arg  == 'V') uplo  = 2;
  ಯದಿ (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  ಯದಿ (incy == 0)          info = 11;
  ಯದಿ (incx == 0)          info =  8;
  ಯದಿ (lda  < k + 1)       info =  6;
  ಯದಿ (k < 0)              info =  3;
  ಯದಿ (n < 0)              info =  2;
  ಯದಿ (uplo  < 0)          info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order,
	   ಗಣನಾ CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   ಶೂನ್ಯ *VALPHA,
	   ಶೂನ್ಯ  *va, blasint lda,
	   ಶೂನ್ಯ  *vx, blasint incx,
	   ಶೂನ್ಯ *VBETA,
	   ಶೂನ್ಯ  *vy, blasint incy){

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
  ಪೂರ್ಣಾಙ್ಕ uplo;
  blasint info;
#ifdef SMPBUG
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 0;
    ಯದಿ (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ಯದಿ (incy == 0)          info = 11;
    ಯದಿ (incx == 0)          info =  8;
    ಯದಿ (lda  < k + 1)       info =  6;
    ಯದಿ (k < 0)              info =  3;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 3;
    ಯದಿ (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    ಯದಿ (incy == 0)          info = 11;
    ಯದಿ (incx == 0)          info =  8;
    ಯದಿ (lda  < k + 1)       info =  6;
    ಯದಿ (k < 0)              info =  3;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  ಯದಿ (n == 0) ಪ್ರತ್ಯಾಗಮ;

  ಯದಿ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ ((alpha_r == ZERO) && (alpha_i == ZERO)) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  ಯದಿ (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ifdef SMPBUG
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (nthreads == 1) {
#endif

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMPBUG
  } ಅನ್ಯಥಾ {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n,  n * k);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
