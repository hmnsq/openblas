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

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ifdef XDOUBLE
#define ERROR_NAME "QSBMV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DSBMV "
#অন্যথা
#define ERROR_NAME "SSBMV "
#endif

স্থির  পূর্ণাঙ্ক (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
  qsbmv_U, qsbmv_L,
#নোচেদ্যদি defined(DOUBLE)
  dsbmv_U, dsbmv_L,
#অন্যথা
  ssbmv_U, ssbmv_L,
#endif
};

#ifdef SMPTEST
স্থির  পূর্ণাঙ্ক (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
  qsbmv_thread_U, qsbmv_thread_L,
#নোচেদ্যদি defined(DOUBLE)
  dsbmv_thread_U, dsbmv_thread_L,
#অন্যথা
  ssbmv_thread_U, ssbmv_thread_L,
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
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  পূর্ণাঙ্ক uplo;
  FLOAT *buffer;
#ifdef SMPTEST
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

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
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  পূর্ণাঙ্ক uplo;
  blasint info;
#ifdef SMPTEST
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
    যদি (Uplo == CblasUpper)         uplo  = 1;
    যদি (Uplo == CblasLower)         uplo  = 0;

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

  যদি (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), রিক্ত, 0, রিক্ত, 0);

  যদি (alpha == ZERO) প্রত্যাগম;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) x -= (n - 1) * incx;
  যদি (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#ifdef SMPTEST
  nthreads = অঙ্ক_cpu_avail(2);

  যদি (nthreads == 1) {
#endif

  (sbmv[uplo])(n, k, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMPTEST
  } অন্যথা {

    (sbmv_thread[uplo])(n, k, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  প্রত্যাগম;
}
