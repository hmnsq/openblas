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
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् SMP
#यदि_परिभाषितम् __64BIT__
#परिभाषय SMPTEST 1
#यदि_अन्तः
#यदि_अन्तः

#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGER  "
#नोचेद्यदि defined DOUBLE
#परिभाषय ERROR_NAME "DGER  "
#अन्यथा
#परिभाषय ERROR_NAME "SGER  "
#यदि_अन्तः

#परिभाषय GER		GERU_K

#यदि   defined XDOUBLE
#परिभाषय GER_THREAD	qger_thread
#नोचेद्यदि defined DOUBLE
#परिभाषय GER_THREAD	dger_thread
#अन्यथा
#परिभाषय GER_THREAD	sger_thread
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(blasint *M, blasint *N, FLOAT *Alpha,
	  FLOAT *x, blasint *INCX,
	  FLOAT *y, blasint *INCY,
	  FLOAT *a, blasint *LDA){

  blasint    m     = *M;
  blasint    n     = *N;
  FLOAT  alpha = *Alpha;
  blasint    incx  = *INCX;
  blasint    incy  = *INCY;
  blasint    lda   = *LDA;
  FLOAT *बफर;
#यदि_परिभाषितम् SMPTEST
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  यदि (lda < MAX(1,m)) info = 9;
  यदि (incy == 0)      info = 7;
  यदि (incx == 0)      info = 5;
  यदि (n < 0)          info = 2;
  यदि (m < 0)          info = 1;

  यदि (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a, blasint lda) {

  FLOAT *बफर;
  blasint info, t;
#यदि_परिभाषितम् SMPTEST
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    info = -1;

    यदि (lda < MAX(1,m)) info = 9;
    यदि (incy == 0)      info = 7;
    यदि (incx == 0)      info = 5;
    यदि (n < 0)          info = 2;
    यदि (m < 0)          info = 1;
  }

  यदि (क्रमः == CblasRowMajor) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    बफर = x;
    x = y;
    y = बफर;

    यदि (lda < MAX(1,m)) info = 9;
    यदि (incy == 0)      info = 7;
    यदि (incx == 0)      info = 5;
    यदि (n < 0)          info = 2;
    यदि (m < 0)          info = 1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  /*     Quick return if possible. */
  यदि (m == 0 || n == 0) प्रत्यागम;
  यदि (alpha == 0.) प्रत्यागम;

  यदि (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_देहली) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, रिक्त);
    प्रत्यागम;
  }  

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incy < 0) y -= (n - 1) * incy;
  यदि (incx < 0) x -= (m - 1) * incx;

  राशिः_ALLOC(m, FLOAT, बफर);

#यदि_परिभाषितम् SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  यदि(1L * m * n > 20480L * GEMM_MULTITHREAD_देहली)
    nthreads = अङ्कः_cpu_avail(2);
  अन्यथा यदि(1L * m * n > 2048L * GEMM_MULTITHREAD_देहली)
	  nthreads = MIN(3,अङ्कः_cpu_avail(2));
  अन्यथा
    nthreads = 1;

  यदि (nthreads == 1) {
#यदि_अन्तः

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, बफर);

#यदि_परिभाषितम् SMPTEST
  } अन्यथा {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, बफर, nthreads);

  }
#यदि_अन्तः

  राशिः_FREE(बफर);
  FUNCTION_PROFILE_अन्तः(1, m * n + m + n, 2 * m * n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
