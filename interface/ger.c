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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGER  "
#elif defined DOUBLE
#define ERROR_NAME "DGER  "
#иначе
#define ERROR_NAME "SGER  "
#endif

#define GER		GERU_K

#если   defined XDOUBLE
#define GER_THREAD	qger_thread
#elif defined DOUBLE
#define GER_THREAD	dger_thread
#иначе
#define GER_THREAD	sger_thread
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(blasint *M, blasint *N, FLOAT *Alpha,
	  FLOAT *x, blasint *INCX,
	  FLOAT *y, blasint *INCY,
	  FLOAT *a, blasint *LDA){

  blasint    m     = *M;
  blasint    n     = *N;
  FLOAT  alpha = *Alpha;
  blasint    incx  = *INCX;
  blasint    incy  = *INCY;
  blasint    lda   = *LDA;
  FLOAT *buffer;
#ifdef SMPTEST
  целое nthreads;
#endif

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  если (lda < MAX(1,m)) info = 9;
  если (incy == 0)      info = 7;
  если (incx == 0)      info = 5;
  если (n < 0)          info = 2;
  если (m < 0)          info = 1;

  если (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#иначе

void CNAME(enum CBLAS_ORDER order,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a, blasint lda) {

  FLOAT *buffer;
  blasint info, t;
#ifdef SMPTEST
  целое nthreads;
#endif

  PRINT_DEBUG_CNAME;

  info  =  0;

  если (order == CblasColMajor) {
    info = -1;

    если (lda < MAX(1,m)) info = 9;
    если (incy == 0)      info = 7;
    если (incx == 0)      info = 5;
    если (n < 0)          info = 2;
    если (m < 0)          info = 1;
  }

  если (order == CblasRowMajor) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    buffer = x;
    x = y;
    y = buffer;

    если (lda < MAX(1,m)) info = 9;
    если (incy == 0)      info = 7;
    если (incx == 0)      info = 5;
    если (n < 0)          info = 2;
    если (m < 0)          info = 1;
  }

  если (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#endif

  /*     Quick return if possible. */
  если (m == 0 || n == 0) вернуть;
  если (alpha == 0.) вернуть;

  если (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_THRESHOLD) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, NULL);
    вернуть;
  }  

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  если (incy < 0) y -= (n - 1) * incy;
  если (incx < 0) x -= (m - 1) * incx;

  STACK_ALLOC(m, FLOAT, buffer);

#ifdef SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  если(1L * m * n > 20480L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = num_cpu_avail(2);
  иначе если(1L * m * n > 2048L * GEMM_MULTITHREAD_THRESHOLD)
	  nthreads = MIN(3,num_cpu_avail(2));
  иначе
    nthreads = 1;

  если (nthreads == 1) {
#endif

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, buffer);

#ifdef SMPTEST
  } иначе {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  STACK_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n, 2 * m * n);

  IDEBUG_END;

  вернуть;
}
