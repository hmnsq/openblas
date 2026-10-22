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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ __64BIT__
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMPTEST 1
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGER  "
#ਜਾਂਜੇ defined DOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGER  "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGER  "
#ਜੇ_ਅੰਤ

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GER		GERU_K

#ਜੇ   defined XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GER_THREAD	qger_thread
#ਜਾਂਜੇ defined DOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GER_THREAD	dger_thread
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GER_THREAD	sger_thread
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(blasint *M, blasint *N, FLOAT *Alpha,
	  FLOAT *x, blasint *INCX,
	  FLOAT *y, blasint *INCY,
	  FLOAT *a, blasint *LDA){

  blasint    m     = *M;
  blasint    n     = *N;
  FLOAT  alpha = *Alpha;
  blasint    incx  = *INCX;
  blasint    incy  = *INCY;
  blasint    lda   = *LDA;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPTEST
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  ਜੇ (lda < MAX(1,m)) info = 9;
  ਜੇ (incy == 0)      info = 7;
  ਜੇ (incx == 0)      info = 5;
  ਜੇ (n < 0)          info = 2;
  ਜੇ (m < 0)          info = 1;

  ਜੇ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a, blasint lda) {

  FLOAT *ਬਫਰ;
  blasint info, t;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPTEST
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    info = -1;

    ਜੇ (lda < MAX(1,m)) info = 9;
    ਜੇ (incy == 0)      info = 7;
    ਜੇ (incx == 0)      info = 5;
    ਜੇ (n < 0)          info = 2;
    ਜੇ (m < 0)          info = 1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    ਬਫਰ = x;
    x = y;
    y = ਬਫਰ;

    ਜੇ (lda < MAX(1,m)) info = 9;
    ਜੇ (incy == 0)      info = 7;
    ਜੇ (incx == 0)      info = 5;
    ਜੇ (n < 0)          info = 2;
    ਜੇ (m < 0)          info = 1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  /*     Quick return if possible. */
  ਜੇ (m == 0 || n == 0) ਵਾਪਸ;
  ਜੇ (alpha == 0.) ਵਾਪਸ;

  ਜੇ (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_ਹੱਦ) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, ਨਿਰਾ);
    ਵਾਪਸ;
  }  

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incy < 0) y -= (n - 1) * incy;
  ਜੇ (incx < 0) x -= (m - 1) * incx;

  ਢੇਰ_ALLOC(m, FLOAT, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  ਜੇ(1L * m * n > 20480L * GEMM_MULTITHREAD_ਹੱਦ)
    nthreads = ਅੰਕ_cpu_avail(2);
  ਨਹੀਂਤਾਂ ਜੇ(1L * m * n > 2048L * GEMM_MULTITHREAD_ਹੱਦ)
	  nthreads = MIN(3,ਅੰਕ_cpu_avail(2));
  ਨਹੀਂਤਾਂ
    nthreads = 1;

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMPTEST
  } ਨਹੀਂਤਾਂ {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  ਢੇਰ_FREE(ਬਫਰ);
  FUNCTION_PROFILE_ਸਿਰਾ(1, m * n + m + n, 2 * m * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
