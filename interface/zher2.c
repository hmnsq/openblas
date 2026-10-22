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
#define ERROR_NAME "XHER2 "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZHER2 "
#অন্যথা
#define ERROR_NAME "CHER2 "
#endif

স্থির পূর্ণাঙ্ক (*her2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xher2_U, xher2_L, xher2_V, xher2_M,
#নোচেদ্যদি defined(DOUBLE)
  zher2_U, zher2_L, zher2_V, zher2_M,
#অন্যথা
  cher2_U, cher2_L, cher2_V, cher2_M,
#endif
};

#ifdef SMP
স্থির পূর্ণাঙ্ক (*her2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
  xher2_thread_U, xher2_thread_L, xher2_thread_V, xher2_thread_M,
#নোচেদ্যদি defined(DOUBLE)
  zher2_thread_U, zher2_thread_L, zher2_thread_V, zher2_thread_M,
#অন্যথা
  cher2_thread_U, cher2_thread_L, cher2_thread_V, cher2_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(বর্ণ *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a, blasint *LDA){

  বর্ণ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  পূর্ণাঙ্ক  uplo;
  FLOAT *buffer;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  যদি (lda  < MAX(1, n))   info =  9;
  যদি (incy == 0)          info =  7;
  যদি (incx == 0)          info =  5;
  যদি (n < 0)              info =  2;
  যদি (uplo  < 0)          info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo, blasint n, শূন্য *VALPHA, শূন্য *vx, blasint incx, শূন্য *vy, blasint incy, শূন্য *va, blasint lda) {

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT *buffer;
  পূর্ণাঙ্ক uplo;
  blasint info;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  যদি (order == CblasColMajor) {

    যদি (Uplo == CblasUpper) uplo  = 0;
    যদি (Uplo == CblasLower) uplo  = 1;

    info = -1;

    যদি (lda  < MAX(1, n))   info =  9;
    যদি (incy == 0)          info =  7;
    যদি (incx == 0)          info =  5;
    যদি (n < 0)              info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (order == CblasRowMajor) {

    যদি (Uplo == CblasUpper) uplo  = 3;
    যদি (Uplo == CblasLower) uplo  = 2;

    info = -1;

    যদি (lda  < MAX(1, n))   info =  9;
    যদি (incx == 0)          info =  7;
    যদি (incy == 0)          info =  5;
    যদি (n < 0)              info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#endif

  যদি (n == 0) প্রত্যাগম;

  যদি ((alpha_r == ZERO) && (alpha_i == ZERO)) প্রত্যাগম;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) x -= (n - 1) * incx * 2;
  যদি (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#ifdef SMP
  যদি (n < 100)
	nthreads = 1;
  অন্যথা
    nthreads = অঙ্ক_cpu_avail(2);

  যদি (nthreads == 1) {
#endif

    (her2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);

#ifdef SMP
  } অন্যথা {

    (her2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  প্রত্যাগম;
}
