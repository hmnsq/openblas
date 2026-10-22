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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_ٿيل __64BIT__
#وصف_ڪريو SMPTEST 1
#جيڪڏهن_پورو
#جيڪڏهن_پورو


#جيڪڏهن_وصف_ٿيل XDOUBLE
#جيڪڏهن_وصف_نه_ٿيل CONJ
#وصف_ڪريو ERROR_NAME "XGERU  "
#وإلا
#وصف_ڪريو ERROR_NAME "XGERC  "
#جيڪڏهن_پورو
#وإلاإذا defined DOUBLE
#جيڪڏهن_وصف_نه_ٿيل CONJ
#وصف_ڪريو ERROR_NAME "ZGERU  "
#وإلا
#وصف_ڪريو ERROR_NAME "ZGERC  "
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_نه_ٿيل CONJ
#وصف_ڪريو ERROR_NAME "CGERU  "
#وإلا
#وصف_ڪريو ERROR_NAME "CGERC "
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#إذا   defined XDOUBLE
#جيڪڏهن_وصف_نه_ٿيل CONJ
#وصف_ڪريو GER		GERU_K
#وصف_ڪريو GER_THREAD	xger_thread_U
#وإلا
#وصف_ڪريو GER		GERC_K
#وصف_ڪريو GER_THREAD	xger_thread_C
#وصف_ڪريو GERV		GERV_K
#وصف_ڪريو GERV_THREAD	xger_thread_V
#جيڪڏهن_پورو
#وإلاإذا defined DOUBLE
#جيڪڏهن_وصف_نه_ٿيل CONJ
#وصف_ڪريو GER		GERU_K
#وصف_ڪريو GER_THREAD	zger_thread_U
#وإلا
#وصف_ڪريو GER		GERC_K
#وصف_ڪريو GER_THREAD	zger_thread_C
#وصف_ڪريو GERV		GERV_K
#وصف_ڪريو GERV_THREAD	zger_thread_V
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_نه_ٿيل CONJ
#وصف_ڪريو GER		GERU_K
#وصف_ڪريو GER_THREAD	cger_thread_U
#وإلا
#وصف_ڪريو GER		GERC_K
#وصف_ڪريو GER_THREAD	cger_thread_C
#وصف_ڪريو GERV		GERV_K
#وصف_ڪريو GERV_THREAD	cger_thread_V
#جيڪڏهن_پورو
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

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
#جيڪڏهن_وصف_ٿيل SMPTEST
  صحيح nthreads;
#جيڪڏهن_پورو

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

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب,
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
#جيڪڏهن_وصف_ٿيل SMPTEST
  صحيح nthreads;
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    info = -1;

    إذا (lda < MAX(1,m)) info = 9;
    إذا (incy == 0)      info = 7;
    إذا (incx == 0)      info = 5;
    إذا (n < 0)          info = 2;
    إذا (m < 0)          info = 1;
  }

  إذا (ترتيب == Cblas_صف_Major) {
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

#جيڪڏهن_پورو

  /*     Quick return if possible. */
  إذا (m == 0 || n == 0) أرجع;

  إذا ((alpha_r == 0.) && (alpha_i == 0.)) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incy < 0) y -= (n - 1) * incy * 2;
  إذا (incx < 0) x -= (m - 1) * incx * 2;

  ڍير_ALLOC(2 * m, FLOAT, بفر);

#جيڪڏهن_وصف_ٿيل SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  إذا(1L * m * n > 36L * حجم(FLOAT) * حجم(FLOAT) * GEMM_MULTITHREAD_حد)
    nthreads = انگ_cpu_avail(2);
  وإلا
    nthreads = 1;

  إذا (nthreads == 1) {
#جيڪڏهن_پورو

#إذا !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
#وإلا
  إذا (ترتيب == Cblas_کانو_Major) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
  } وإلا {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);
  }
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMPTEST

  } وإلا {

#إذا !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
#وإلا
      إذا (ترتيب == Cblas_کانو_Major) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
      } وإلا {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, بفر, nthreads);
      }
#جيڪڏهن_پورو

  }
#جيڪڏهن_پورو

  ڍير_FREE(بفر);

  FUNCTION_PROFILE_پڇاڙي(4, m * n + m + n, 2 * m * n);

  IDEBUG_پڇاڙي;

  أرجع;

}
