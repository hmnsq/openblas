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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHBMV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZHBMV "
#অন্যথা
#define ERROR_NAME "CHBMV "
#endif

স্থির  পূর্ণাঙ্ক (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#নোচেদ্যদি defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#অন্যথা
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#endif
};

#ifdef SMPBUG
স্থির  পূর্ণাঙ্ক (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#নোচেদ্যদি defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#অন্যথা
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(বর্ণ *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  বর্ণ uplo_arg = *UPLO;
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
  পূর্ণাঙ্ক uplo;
  FLOAT *buffer;
#ifdef SMPBUG
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;
  যদি (uplo_arg  == 'V') uplo  = 2;
  যদি (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  যদি (incy == 0)          info = 11;
  যদি (incx == 0)          info =  8;
  যদি (lda  < k + 1)       info =  6;
  যদি (k < 0)              info =  3;
  যদি (n < 0)              info =  2;
  যদি (uplo  < 0)          info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order,
	   গণনা CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   শূন্য *VALPHA,
	   শূন্য  *va, blasint lda,
	   শূন্য  *vx, blasint incx,
	   শূন্য *VBETA,
	   শূন্য  *vy, blasint incy){

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
  পূর্ণাঙ্ক uplo;
  blasint info;
#ifdef SMPBUG
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (Uplo == CblasUpper)         uplo  = 0;
    যদি (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    যদি (incy == 0)          info = 11;
    যদি (incx == 0)          info =  8;
    যদি (lda  < k + 1)       info =  6;
    যদি (k < 0)              info =  3;
    যদি (n < 0)              info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (order == CblasRowMajor) {
    যদি (Uplo == CblasUpper)         uplo  = 3;
    যদি (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    যদি (incy == 0)          info = 11;
    যদি (incx == 0)          info =  8;
    যদি (lda  < k + 1)       info =  6;
    যদি (k < 0)              info =  3;
    যদি (n < 0)              info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#endif

  যদি (n == 0) প্রত্যাগম;

  যদি ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), রিক্ত, 0, রিক্ত, 0);

  যদি ((alpha_r == ZERO) && (alpha_i == ZERO)) প্রত্যাগম;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  যদি (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#ifdef SMPBUG
  nthreads = অঙ্ক_cpu_avail(2);

  যদি (nthreads == 1) {
#endif

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMPBUG
  } অন্যথা {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n,  n * k);

  IDEBUG_END;

  প্রত্যাগম;
}
