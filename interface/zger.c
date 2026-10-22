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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif


#ifdef XDOUBLE
#ifndef CONJ
#define ERROR_NAME "XGERU  "
#अञ्ञथा
#define ERROR_NAME "XGERC  "
#endif
#नोचेद्यदि defined DOUBLE
#ifndef CONJ
#define ERROR_NAME "ZGERU  "
#अञ्ञथा
#define ERROR_NAME "ZGERC  "
#endif
#अञ्ञथा
#ifndef CONJ
#define ERROR_NAME "CGERU  "
#अञ्ञथा
#define ERROR_NAME "CGERC "
#endif
#endif

#सचे   defined XDOUBLE
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	xger_thread_U
#अञ्ञथा
#define GER		GERC_K
#define GER_THREAD	xger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	xger_thread_V
#endif
#नोचेद्यदि defined DOUBLE
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	zger_thread_U
#अञ्ञथा
#define GER		GERC_K
#define GER_THREAD	zger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	zger_thread_V
#endif
#अञ्ञथा
#ifndef CONJ
#define GER		GERU_K
#define GER_THREAD	cger_thread_U
#अञ्ञथा
#define GER		GERC_K
#define GER_THREAD	cger_thread_C
#define GERV		GERV_K
#define GERV_THREAD	cger_thread_V
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
  पूर्णाङ्क nthreads;
#endif

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  सचे (lda < MAX(1,m)) info = 9;
  सचे (incy == 0)      info = 7;
  सचे (incx == 0)      info = 5;
  सचे (n < 0)          info = 2;
  सचे (m < 0)          info = 1;

  सचे (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   blasint m, blasint n,
	   शून्य *VAlpha,
	   शून्य  *vx, blasint incx,
	   शून्य  *vy, blasint incy,
	   शून्य  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *buffer;
  blasint info, t;
#ifdef SMPTEST
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  info  =  0;

  सचे (order == CblasColMajor) {
    info = -1;

    सचे (lda < MAX(1,m)) info = 9;
    सचे (incy == 0)      info = 7;
    सचे (incx == 0)      info = 5;
    सचे (n < 0)          info = 2;
    सचे (m < 0)          info = 1;
  }

  सचे (order == CblasRowMajor) {
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

    सचे (lda < MAX(1,m)) info = 9;
    सचे (incy == 0)      info = 7;
    सचे (incx == 0)      info = 5;
    सचे (n < 0)          info = 2;
    सचे (m < 0)          info = 1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  /*     Quick return if possible. */
  सचे (m == 0 || n == 0) पटिनिवत्त;

  सचे ((alpha_r == 0.) && (alpha_i == 0.)) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incy < 0) y -= (n - 1) * incy * 2;
  सचे (incx < 0) x -= (m - 1) * incx * 2;

  राशिः_ALLOC(2 * m, FLOAT, buffer);

#ifdef SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  सचे(1L * m * n > 36L * परिमाण(FLOAT) * परिमाण(FLOAT) * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = अङ्कः_cpu_avail(2);
  अञ्ञथा
    nthreads = 1;

  सचे (nthreads == 1) {
#endif

#सचे !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
#अञ्ञथा
  सचे (order == CblasColMajor) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
  } अञ्ञथा {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);
  }
#endif

#ifdef SMPTEST

  } अञ्ञथा {

#सचे !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
#अञ्ञथा
      सचे (order == CblasColMajor) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
      } अञ्ञथा {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, buffer, nthreads);
      }
#endif

  }
#endif

  राशिः_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n, 2 * m * n);

  IDEBUG_END;

  पटिनिवत्त;

}
