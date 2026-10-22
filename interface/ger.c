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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGER  "
#અથવાજો defined DOUBLE
#define ERROR_NAME "DGER  "
#નહીંતર
#define ERROR_NAME "SGER  "
#endif

#define GER		GERU_K

#જો   defined XDOUBLE
#define GER_THREAD	qger_thread
#અથવાજો defined DOUBLE
#define GER_THREAD	dger_thread
#નહીંતર
#define GER_THREAD	sger_thread
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
  પૂર્ણાંક nthreads;
#endif

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  જો (lda < MAX(1,m)) info = 9;
  જો (incy == 0)      info = 7;
  જો (incx == 0)      info = 5;
  જો (n < 0)          info = 2;
  જો (m < 0)          info = 1;

  જો (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a, blasint lda) {

  FLOAT *buffer;
  blasint info, t;
#ifdef SMPTEST
  પૂર્ણાંક nthreads;
#endif

  PRINT_DEBUG_CNAME;

  info  =  0;

  જો (order == CblasColMajor) {
    info = -1;

    જો (lda < MAX(1,m)) info = 9;
    જો (incy == 0)      info = 7;
    જો (incx == 0)      info = 5;
    જો (n < 0)          info = 2;
    જો (m < 0)          info = 1;
  }

  જો (order == CblasRowMajor) {
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

    જો (lda < MAX(1,m)) info = 9;
    જો (incy == 0)      info = 7;
    જો (incx == 0)      info = 5;
    જો (n < 0)          info = 2;
    જો (m < 0)          info = 1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  /*     Quick return if possible. */
  જો (m == 0 || n == 0) પાછા;
  જો (alpha == 0.) પાછા;

  જો (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_THRESHOLD) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, ખાલી);
    પાછા;
  }  

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incy < 0) y -= (n - 1) * incy;
  જો (incx < 0) x -= (m - 1) * incx;

  ઢગલો_ALLOC(m, FLOAT, buffer);

#ifdef SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  જો(1L * m * n > 20480L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = અંક_cpu_avail(2);
  નહીંતર જો(1L * m * n > 2048L * GEMM_MULTITHREAD_THRESHOLD)
	  nthreads = MIN(3,અંક_cpu_avail(2));
  નહીંતર
    nthreads = 1;

  જો (nthreads == 1) {
#endif

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, buffer);

#ifdef SMPTEST
  } નહીંતર {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  ઢગલો_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n, 2 * m * n);

  IDEBUG_END;

  પાછા;
}
