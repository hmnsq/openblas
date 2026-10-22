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

#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ SMP
#اگر_تعریف_شدہ __64BIT__
#تعریف_کرو SMPTEST 1
#اگر_ختم
#اگر_ختم


#اگر_تعریف_شدہ XDOUBLE
#اگر_غیر_تعریف_شدہ CONJ
#تعریف_کرو ERROR_NAME "XGERU  "
#ورنہ
#تعریف_کرو ERROR_NAME "XGERC  "
#اگر_ختم
#وگرنہ_اگر defined DOUBLE
#اگر_غیر_تعریف_شدہ CONJ
#تعریف_کرو ERROR_NAME "ZGERU  "
#ورنہ
#تعریف_کرو ERROR_NAME "ZGERC  "
#اگر_ختم
#ورنہ
#اگر_غیر_تعریف_شدہ CONJ
#تعریف_کرو ERROR_NAME "CGERU  "
#ورنہ
#تعریف_کرو ERROR_NAME "CGERC "
#اگر_ختم
#اگر_ختم

#اگر   defined XDOUBLE
#اگر_غیر_تعریف_شدہ CONJ
#تعریف_کرو GER		GERU_K
#تعریف_کرو GER_THREAD	xger_thread_U
#ورنہ
#تعریف_کرو GER		GERC_K
#تعریف_کرو GER_THREAD	xger_thread_C
#تعریف_کرو GERV		GERV_K
#تعریف_کرو GERV_THREAD	xger_thread_V
#اگر_ختم
#وگرنہ_اگر defined DOUBLE
#اگر_غیر_تعریف_شدہ CONJ
#تعریف_کرو GER		GERU_K
#تعریف_کرو GER_THREAD	zger_thread_U
#ورنہ
#تعریف_کرو GER		GERC_K
#تعریف_کرو GER_THREAD	zger_thread_C
#تعریف_کرو GERV		GERV_K
#تعریف_کرو GERV_THREAD	zger_thread_V
#اگر_ختم
#ورنہ
#اگر_غیر_تعریف_شدہ CONJ
#تعریف_کرو GER		GERU_K
#تعریف_کرو GER_THREAD	cger_thread_U
#ورنہ
#تعریف_کرو GER		GERC_K
#تعریف_کرو GER_THREAD	cger_thread_C
#تعریف_کرو GERV		GERV_K
#تعریف_کرو GERV_THREAD	cger_thread_V
#اگر_ختم
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
  FLOAT *بفر;
#اگر_تعریف_شدہ SMPTEST
  صحیح nthreads;
#اگر_ختم

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  اگر (lda < MAX(1,m)) info = 9;
  اگر (incy == 0)      info = 7;
  اگر (incx == 0)      info = 5;
  اگر (n < 0)          info = 2;
  اگر (m < 0)          info = 1;

  اگر (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب,
	   blasint m, blasint n,
	   خالی *VAlpha,
	   خالی  *vx, blasint incx,
	   خالی  *vy, blasint incy,
	   خالی  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *بفر;
  blasint info, t;
#اگر_تعریف_شدہ SMPTEST
  صحیح nthreads;
#اگر_ختم

  PRINT_DEBUG_CNAME;

  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    info = -1;

    اگر (lda < MAX(1,m)) info = 9;
    اگر (incy == 0)      info = 7;
    اگر (incx == 0)      info = 5;
    اگر (n < 0)          info = 2;
    اگر (m < 0)          info = 1;
  }

  اگر (ترتیب == Cblas_صف_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    بفر = x;
    x = y;
    y = بفر;

    اگر (lda < MAX(1,m)) info = 9;
    اگر (incy == 0)      info = 7;
    اگر (incx == 0)      info = 5;
    اگر (n < 0)          info = 2;
    اگر (m < 0)          info = 1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_ختم

  /*     Quick return if possible. */
  اگر (m == 0 || n == 0) واپس;

  اگر ((alpha_r == 0.) && (alpha_i == 0.)) واپس;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incy < 0) y -= (n - 1) * incy * 2;
  اگر (incx < 0) x -= (m - 1) * incx * 2;

  ڈھیر_ALLOC(2 * m, FLOAT, بفر);

#اگر_تعریف_شدہ SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  اگر(1L * m * n > 36L * حجم(FLOAT) * حجم(FLOAT) * GEMM_MULTITHREAD_حد)
    nthreads = عدد_cpu_avail(2);
  ورنہ
    nthreads = 1;

  اگر (nthreads == 1) {
#اگر_ختم

#اگر !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
#ورنہ
  اگر (ترتیب == Cblas_ستون_Major) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
  } ورنہ {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
  }
#اگر_ختم

#اگر_تعریف_شدہ SMPTEST

  } ورنہ {

#اگر !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
#ورنہ
      اگر (ترتیب == Cblas_ستون_Major) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
      } ورنہ {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
      }
#اگر_ختم

  }
#اگر_ختم

  ڈھیر_FREE(بفر);

  FUNCTION_PROFILE_انتہا(4, m * n + m + n, 2 * m * n);

  IDEBUG_انتہا;

  واپس;

}
