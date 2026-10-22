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
#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_تعریف_شدہ SMP
#جے_تعریف_شدہ __64BIT__
#تعریف_کرو SMPTEST 1
#جے_انت
#جے_انت


#جے_تعریف_شدہ XDOUBLE
#جے_بنا_تعریف CONJ
#تعریف_کرو ERROR_NAME "XGERU  "
#وإلا
#تعریف_کرو ERROR_NAME "XGERC  "
#جے_انت
#وإلاإذا defined DOUBLE
#جے_بنا_تعریف CONJ
#تعریف_کرو ERROR_NAME "ZGERU  "
#وإلا
#تعریف_کرو ERROR_NAME "ZGERC  "
#جے_انت
#وإلا
#جے_بنا_تعریف CONJ
#تعریف_کرو ERROR_NAME "CGERU  "
#وإلا
#تعریف_کرو ERROR_NAME "CGERC "
#جے_انت
#جے_انت

#إذا   defined XDOUBLE
#جے_بنا_تعریف CONJ
#تعریف_کرو GER		GERU_K
#تعریف_کرو GER_THREAD	xger_thread_U
#وإلا
#تعریف_کرو GER		GERC_K
#تعریف_کرو GER_THREAD	xger_thread_C
#تعریف_کرو GERV		GERV_K
#تعریف_کرو GERV_THREAD	xger_thread_V
#جے_انت
#وإلاإذا defined DOUBLE
#جے_بنا_تعریف CONJ
#تعریف_کرو GER		GERU_K
#تعریف_کرو GER_THREAD	zger_thread_U
#وإلا
#تعریف_کرو GER		GERC_K
#تعریف_کرو GER_THREAD	zger_thread_C
#تعریف_کرو GERV		GERV_K
#تعریف_کرو GERV_THREAD	zger_thread_V
#جے_انت
#وإلا
#جے_بنا_تعریف CONJ
#تعریف_کرو GER		GERU_K
#تعریف_کرو GER_THREAD	cger_thread_U
#وإلا
#تعریف_کرو GER		GERC_K
#تعریف_کرو GER_THREAD	cger_thread_C
#تعریف_کرو GERV		GERV_K
#تعریف_کرو GERV_THREAD	cger_thread_V
#جے_انت
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
#جے_تعریف_شدہ SMPTEST
  صحيح nthreads;
#جے_انت

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  إذا (lda < MAX(1,m)) info = 9;
  إذا (incy == 0)      info = 7;
  إذا (incx == 0)      info = 5;
  إذا (n < 0)          info = 2;
  إذا (m < 0)          info = 1;

  إذا (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب,
	   blasint m, blasint n,
	   فراغ *VAlpha,
	   فراغ  *vx, blasint incx,
	   فراغ  *vy, blasint incy,
	   فراغ  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *بفر;
  blasint info, t;
#جے_تعریف_شدہ SMPTEST
  صحيح nthreads;
#جے_انت

  PRINT_DEBUG_CNAME;

  info  =  0;

  إذا (ترتیب == Cblas_خانہ_Major) {
    info = -1;

    إذا (lda < MAX(1,m)) info = 9;
    إذا (incy == 0)      info = 7;
    إذا (incx == 0)      info = 5;
    إذا (n < 0)          info = 2;
    إذا (m < 0)          info = 1;
  }

  إذا (ترتیب == Cblas_صف_Major) {
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

    إذا (lda < MAX(1,m)) info = 9;
    إذا (incy == 0)      info = 7;
    إذا (incx == 0)      info = 5;
    إذا (n < 0)          info = 2;
    إذا (m < 0)          info = 1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جے_انت

  /*     Quick return if possible. */
  إذا (m == 0 || n == 0) أرجع;

  إذا ((alpha_r == 0.) && (alpha_i == 0.)) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incy < 0) y -= (n - 1) * incy * 2;
  إذا (incx < 0) x -= (m - 1) * incx * 2;

  ڈھیر_ALLOC(2 * m, FLOAT, بفر);

#جے_تعریف_شدہ SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  إذا(1L * m * n > 36L * حجم(FLOAT) * حجم(FLOAT) * GEMM_MULTITHREAD_حد)
    nthreads = عدد_cpu_avail(2);
  وإلا
    nthreads = 1;

  إذا (nthreads == 1) {
#جے_انت

#إذا !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
#وإلا
  إذا (ترتیب == Cblas_خانہ_Major) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
  } وإلا {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
  }
#جے_انت

#جے_تعریف_شدہ SMPTEST

  } وإلا {

#إذا !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
#وإلا
      إذا (ترتیب == Cblas_خانہ_Major) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
      } وإلا {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
      }
#جے_انت

  }
#جے_انت

  ڈھیر_FREE(بفر);

  FUNCTION_PROFILE_سرا(4, m * n + m + n, 2 * m * n);

  IDEBUG_سرا;

  أرجع;

}
