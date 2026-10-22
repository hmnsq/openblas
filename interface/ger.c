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
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_شده SMP
#اگر_تعریف_شده __64BIT__
#تعریف_کن SMPTEST 1
#پایان_اگر
#پایان_اگر

#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QGER  "
#وإلاإذا defined DOUBLE
#تعریف_کن ERROR_NAME "DGER  "
#وإلا
#تعریف_کن ERROR_NAME "SGER  "
#پایان_اگر

#تعریف_کن GER		GERU_K

#إذا   defined XDOUBLE
#تعریف_کن GER_THREAD	qger_thread
#وإلاإذا defined DOUBLE
#تعریف_کن GER_THREAD	dger_thread
#وإلا
#تعریف_کن GER_THREAD	sger_thread
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(blasint *M, blasint *N, FLOAT *Alpha,
	  FLOAT *x, blasint *INCX,
	  FLOAT *y, blasint *INCY,
	  FLOAT *a, blasint *LDA){

  blasint    m     = *M;
  blasint    n     = *N;
  FLOAT  alpha = *Alpha;
  blasint    incx  = *INCX;
  blasint    incy  = *INCY;
  blasint    lda   = *LDA;
  FLOAT *میانگیر;
#اگر_تعریف_شده SMPTEST
  صحيح nthreads;
#پایان_اگر

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
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a, blasint lda) {

  FLOAT *میانگیر;
  blasint info, t;
#اگر_تعریف_شده SMPTEST
  صحيح nthreads;
#پایان_اگر

  PRINT_DEBUG_CNAME;

  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    info = -1;

    إذا (lda < MAX(1,m)) info = 9;
    إذا (incy == 0)      info = 7;
    إذا (incx == 0)      info = 5;
    إذا (n < 0)          info = 2;
    إذا (m < 0)          info = 1;
  }

  إذا (ترتیب == Cblas_ردیف_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    میانگیر = x;
    x = y;
    y = میانگیر;

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

#پایان_اگر

  /*     Quick return if possible. */
  إذا (m == 0 || n == 0) أرجع;
  إذا (alpha == 0.) أرجع;

  إذا (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_آستانه) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, عدم);
    أرجع;
  }  

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incy < 0) y -= (n - 1) * incy;
  إذا (incx < 0) x -= (m - 1) * incx;

  پشته_ALLOC(m, FLOAT, میانگیر);

#اگر_تعریف_شده SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  إذا(1L * m * n > 20480L * GEMM_MULTITHREAD_آستانه)
    nthreads = رقم_cpu_avail(2);
  وإلا إذا(1L * m * n > 2048L * GEMM_MULTITHREAD_آستانه)
	  nthreads = MIN(3,رقم_cpu_avail(2));
  وإلا
    nthreads = 1;

  إذا (nthreads == 1) {
#پایان_اگر

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, میانگیر);

#اگر_تعریف_شده SMPTEST
  } وإلا {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, میانگیر, nthreads);

  }
#پایان_اگر

  پشته_FREE(میانگیر);
  FUNCTION_PROFILE_انتها(1, m * n + m + n, 2 * m * n);

  IDEBUG_انتها;

  أرجع;
}
