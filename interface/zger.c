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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif


#ifdef XDOUBLE
#ifndef CONJ
#define ERROR_NAME "XGERU  "
#இல்லையெனில்
#define ERROR_NAME "XGERC  "
#endif
#அல்லதெனில் defined DOUBLE
#ifndef CONJ
#define ERROR_NAME "ZGERU  "
#இல்லையெனில்
#define ERROR_NAME "ZGERC  "
#endif
#இல்லையெனில்
#ifndef CONJ
#define ERROR_NAME "CGERU  "
#இல்லையெனில்
#define ERROR_NAME "CGERC "
#endif
#endif

#எனில்   defined XDOUBLE
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	xger_thread_U
#இல்லையெனில்
#define GER		GERC_K
#define GER_THREAD	xger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	xger_thread_V
#endif
#அல்லதெனில் defined DOUBLE
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	zger_thread_U
#இல்லையெனில்
#define GER		GERC_K
#define GER_THREAD	zger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	zger_thread_V
#endif
#இல்லையெனில்
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	cger_thread_U
#இல்லையெனில்
#define GER		GERC_K
#define GER_THREAD	cger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	cger_thread_V
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
  முழு nthreads;
#endif

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  எனில் (lda < MAX(1,m)) info = 9;
  எனில் (incy == 0)      info = 7;
  எனில் (incx == 0)      info = 5;
  எனில் (n < 0)          info = 2;
  எனில் (m < 0)          info = 1;

  எனில் (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order,
	   blasint m, blasint n,
	   வெறுமை *VAlpha,
	   வெறுமை  *vx, blasint incx,
	   வெறுமை  *vy, blasint incy,
	   வெறுமை  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *buffer;
  blasint info, t;
#ifdef SMPTEST
  முழு nthreads;
#endif

  PRINT_DEBUG_CNAME;

  info  =  0;

  எனில் (order == CblasColMajor) {
    info = -1;

    எனில் (lda < MAX(1,m)) info = 9;
    எனில் (incy == 0)      info = 7;
    எனில் (incx == 0)      info = 5;
    எனில் (n < 0)          info = 2;
    எனில் (m < 0)          info = 1;
  }

  எனில் (order == CblasRowMajor) {
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

    எனில் (lda < MAX(1,m)) info = 9;
    எனில் (incy == 0)      info = 7;
    எனில் (incx == 0)      info = 5;
    எனில் (n < 0)          info = 2;
    எனில் (m < 0)          info = 1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#endif

  /*     Quick return if possible. */
  எனில் (m == 0 || n == 0) திரும்பு;

  எனில் ((alpha_r == 0.) && (alpha_i == 0.)) திரும்பு;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  எனில் (incy < 0) y -= (n - 1) * incy * 2;
  எனில் (incx < 0) x -= (m - 1) * incx * 2;

  அடுக்கு_ALLOC(2 * m, FLOAT, buffer);

#ifdef SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  எனில்(1L * m * n > 36L * அளவு(FLOAT) * அளவு(FLOAT) * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = இலக்கம்_cpu_avail(2);
  இல்லையெனில்
    nthreads = 1;

  எனில் (nthreads == 1) {
#endif

#எனில் !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
#இல்லையெனில்
  எனில் (order == CblasColMajor) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
  } இல்லையெனில் {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
  }
#endif

#ifdef SMPTEST

  } இல்லையெனில் {

#எனில் !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
#இல்லையெனில்
      எனில் (order == CblasColMajor) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
      } இல்லையெனில் {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
      }
#endif

  }
#endif

  அடுக்கு_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n, 2 * m * n);

  IDEBUG_END;

  திரும்பு;

}
