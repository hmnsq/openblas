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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_تعریف_شوی SMP
#که_تعریف_شوی __64BIT__
#تعریف_کړه SMPTEST 1
#د_که_پای
#د_که_پای

#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QGER  "
#وإلاإذا defined DOUBLE
#تعریف_کړه ERROR_NAME "DGER  "
#وإلا
#تعریف_کړه ERROR_NAME "SGER  "
#د_که_پای

#تعریف_کړه GER		GERU_K

#إذا   defined XDOUBLE
#تعریف_کړه GER_THREAD	qger_thread
#وإلاإذا defined DOUBLE
#تعریف_کړه GER_THREAD	dger_thread
#وإلا
#تعریف_کړه GER_THREAD	sger_thread
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

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
  FLOAT *بفر;
#که_تعریف_شوی SMPTEST
  صحيح nthreads;
#د_که_پای

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

  FLOAT *بفر;
  blasint info, t;
#که_تعریف_شوی SMPTEST
  صحيح nthreads;
#د_که_پای

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

  إذا (ترتیب == Cblas_لیکه_Major) {
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

#د_که_پای

  /*     Quick return if possible. */
  إذا (m == 0 || n == 0) أرجع;
  إذا (alpha == 0.) أرجع;

  إذا (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_پوله) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, عدم);
    أرجع;
  }  

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incy < 0) y -= (n - 1) * incy;
  إذا (incx < 0) x -= (m - 1) * incx;

  ډېری_ALLOC(m, FLOAT, بفر);

#که_تعریف_شوی SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  إذا(1L * m * n > 20480L * GEMM_MULTITHREAD_پوله)
    nthreads = عدد_cpu_avail(2);
  وإلا إذا(1L * m * n > 2048L * GEMM_MULTITHREAD_پوله)
	  nthreads = MIN(3,عدد_cpu_avail(2));
  وإلا
    nthreads = 1;

  إذا (nthreads == 1) {
#د_که_پای

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, بفر);

#که_تعریف_شوی SMPTEST
  } وإلا {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, بفر, nthreads);

  }
#د_که_پای

  ډېری_FREE(بفر);
  FUNCTION_PROFILE_پای_ټکی(1, m * n + m + n, 2 * m * n);

  IDEBUG_پای_ټکی;

  أرجع;
}
