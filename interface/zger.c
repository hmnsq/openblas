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
#ifndef CONJ
#define ERROR_NAME "XGERU  "
#không_thì
#define ERROR_NAME "XGERC  "
#endif
#elif defined DOUBLE
#ifndef CONJ
#define ERROR_NAME "ZGERU  "
#không_thì
#define ERROR_NAME "ZGERC  "
#endif
#không_thì
#ifndef CONJ
#define ERROR_NAME "CGERU  "
#không_thì
#define ERROR_NAME "CGERC "
#endif
#endif

#nếu   defined XDOUBLE
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	xger_thread_U
#không_thì
#define GER		GERC_K
#define GER_THREAD	xger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	xger_thread_V
#endif
#elif defined DOUBLE
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	zger_thread_U
#không_thì
#define GER		GERC_K
#define GER_THREAD	zger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	zger_thread_V
#endif
#không_thì
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	cger_thread_U
#không_thì
#define GER		GERC_K
#define GER_THREAD	cger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	cger_thread_V
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(blasint *M, blasint *N, FLOAT *Alpha,
	  FLOAT *x, blasint *INCX,
	  FLOAT *y, blasint *INCY,
	  FLOAT *a, blasint *LDA){

  blasint    m     = *M;
  blasint    n     = *N;
  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];
  blasint    incx  = *INCX;
  blasint    incy  = *INCY;
  blasint    lda   = *LDA;
  FLOAT *buffer;
#ifdef SMPTEST
  số_nguyên nthreads;
#endif

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  nếu (lda < MAX(1,m)) info = 9;
  nếu (incy == 0)      info = 7;
  nếu (incx == 0)      info = 5;
  nếu (n < 0)          info = 2;
  nếu (m < 0)          info = 1;

  nếu (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#không_thì

void CNAME(enum CBLAS_ORDER order,
	   blasint m, blasint n,
	   void *VAlpha,
	   void  *vx, blasint incx,
	   void  *vy, blasint incy,
	   void  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *buffer;
  blasint info, t;
#ifdef SMPTEST
  số_nguyên nthreads;
#endif

  PRINT_DEBUG_CNAME;

  info  =  0;

  nếu (order == CblasColMajor) {
    info = -1;

    nếu (lda < MAX(1,m)) info = 9;
    nếu (incy == 0)      info = 7;
    nếu (incx == 0)      info = 5;
    nếu (n < 0)          info = 2;
    nếu (m < 0)          info = 1;
  }

  nếu (order == CblasRowMajor) {
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

    nếu (lda < MAX(1,m)) info = 9;
    nếu (incy == 0)      info = 7;
    nếu (incx == 0)      info = 5;
    nếu (n < 0)          info = 2;
    nếu (m < 0)          info = 1;
  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#endif

  /*     Quick return if possible. */
  nếu (m == 0 || n == 0) trả_về;

  nếu ((alpha_r == 0.) && (alpha_i == 0.)) trả_về;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  nếu (incy < 0) y -= (n - 1) * incy * 2;
  nếu (incx < 0) x -= (m - 1) * incx * 2;

  STACK_ALLOC(2 * m, FLOAT, buffer);

#ifdef SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  nếu(1L * m * n > 36L * sizeof(FLOAT) * sizeof(FLOAT) * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = num_cpu_avail(2);
  không_thì
    nthreads = 1;

  nếu (nthreads == 1) {
#endif

#nếu !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
#không_thì
  nếu (order == CblasColMajor) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
  } không_thì {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
  }
#endif

#ifdef SMPTEST

  } không_thì {

#nếu !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
#không_thì
      nếu (order == CblasColMajor) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
      } không_thì {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
      }
#endif

  }
#endif

  STACK_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n, 2 * m * n);

  IDEBUG_END;

  trả_về;

}
