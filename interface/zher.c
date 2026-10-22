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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHER  "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHER  "
#ܐܚܪܢܐ
#define ERROR_NAME "CHER  "
#endif

ساكن صحيح (*her[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xher_U, xher_L, xher_V, xher_M,
#وإلاإذا defined(DOUBLE)
  zher_U, zher_L, zher_V, zher_M,
#ܐܚܪܢܐ
  cher_U, cher_L, cher_V, cher_M,
#endif
};

#ifdef SMP
ساكن صحيح (*her_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#ifdef XDOUBLE
  xher_thread_U, xher_thread_L, xher_thread_V, xher_thread_M,
#وإلاإذا defined(DOUBLE)
  zher_thread_U, zher_thread_L, zher_thread_V, zher_thread_M,
#ܐܚܪܢܐ
  cher_thread_U, cher_thread_L, cher_thread_V, cher_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	  FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  صحيح uplo;
  FLOAT *buffer;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ܐܢ (uplo_arg  == 'U') uplo  = 0;
  ܐܢ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ܐܢ (lda  < MAX(1, n))   info =  7;
  ܐܢ (incx == 0)          info =  5;
  ܐܢ (n < 0)              info =  2;
  ܐܢ (uplo  < 0)          info =  1;

  ܐܢ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ܐܚܪܢܐ

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo, blasint n, FLOAT alpha, فراغ *vx, blasint incx, فراغ *va, blasint lda) {

  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;

  FLOAT *buffer;
  صحيح uplo;
  blasint info;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ܐܢ (order == CblasColMajor) {

    ܐܢ (Uplo == CblasUpper) uplo  = 0;
    ܐܢ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ܐܢ (lda  < MAX(1, n))   info =  7;
    ܐܢ (incx == 0)          info =  5;
    ܐܢ (n < 0)              info =  2;
    ܐܢ (uplo  < 0)          info =  1;

  }

  ܐܢ (order == CblasRowMajor) {

    ܐܢ (Uplo == CblasUpper) uplo  = 3;
    ܐܢ (Uplo == CblasLower) uplo  = 2;

    info = -1;

    ܐܢ (lda  < MAX(1, n))   info =  7;
    ܐܢ (incx == 0)          info =  5;
    ܐܢ (n < 0)              info =  2;
    ܐܢ (uplo  < 0)          info =  1;
  }

  ܐܢ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#endif

  ܐܢ (n == 0) ܦܢܐ;

  ܐܢ (alpha == ZERO) ܦܢܐ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ܐܢ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ܐܢ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ifdef SMP
  ܐܢ (n < 100)
    nthreads = 1;
  ܐܚܪܢܐ
    nthreads = num_cpu_avail(2);

  ܐܢ (nthreads == 1) {
#endif

  (her[uplo])(n, alpha, x, incx, a, lda, buffer);

#ifdef SMP
  } ܐܚܪܢܐ {

    (her_thread[uplo])(n, alpha, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ܦܢܐ;
}
