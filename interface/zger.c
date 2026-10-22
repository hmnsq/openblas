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
#यदि_अपरिभाषितम् CONJ
#परिभाषय ERROR_NAME "XGERU  "
#अन्यथा
#परिभाषय ERROR_NAME "XGERC  "
#यदि_अन्तः
#नोचेद्यदि defined DOUBLE
#यदि_अपरिभाषितम् CONJ
#परिभाषय ERROR_NAME "ZGERU  "
#अन्यथा
#परिभाषय ERROR_NAME "ZGERC  "
#यदि_अन्तः
#अन्यथा
#यदि_अपरिभाषितम् CONJ
#परिभाषय ERROR_NAME "CGERU  "
#अन्यथा
#परिभाषय ERROR_NAME "CGERC "
#यदि_अन्तः
#यदि_अन्तः

#यदि   defined XDOUBLE
#यदि_अपरिभाषितम् CONJ
#परिभाषय GER		GERU_K
#परिभाषय GER_THREAD	xger_thread_U
#अन्यथा
#परिभाषय GER		GERC_K
#परिभाषय GER_THREAD	xger_thread_C
#परिभाषय GERV		GERV_K
#परिभाषय GERV_THREAD	xger_thread_V
#यदि_अन्तः
#नोचेद्यदि defined DOUBLE
#यदि_अपरिभाषितम् CONJ
#परिभाषय GER		GERU_K
#परिभाषय GER_THREAD	zger_thread_U
#अन्यथा
#परिभाषय GER		GERC_K
#परिभाषय GER_THREAD	zger_thread_C
#परिभाषय GERV		GERV_K
#परिभाषय GERV_THREAD	zger_thread_V
#यदि_अन्तः
#अन्यथा
#यदि_अपरिभाषितम् CONJ
#परिभाषय GER		GERU_K
#परिभाषय GER_THREAD	cger_thread_U
#अन्यथा
#परिभाषय GER		GERC_K
#परिभाषय GER_THREAD	cger_thread_C
#परिभाषय GERV		GERV_K
#परिभाषय GERV_THREAD	cger_thread_V
#यदि_अन्तः
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

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

  यदि ((alpha_r == 0.) && (alpha_i == 0.)) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incy < 0) y -= (n - 1) * incy * 2;
  यदि (incx < 0) x -= (m - 1) * incx * 2;

  राशिः_ALLOC(2 * m, FLOAT, बफर);

#यदि_परिभाषितम् SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  यदि(1L * m * n > 36L * परिमाण(FLOAT) * परिमाण(FLOAT) * GEMM_MULTITHREAD_देहली)
    nthreads = अङ्कः_cpu_avail(2);
  अन्यथा
    nthreads = 1;

  यदि (nthreads == 1) {
#यदि_अन्तः

#यदि !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, बफर);
#अन्यथा
  यदि (क्रमः == CblasColMajor) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, बफर);
  } अन्यथा {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, बफर);
  }
#यदि_अन्तः

#यदि_परिभाषितम् SMPTEST

  } अन्यथा {

#यदि !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, बफर, nthreads);
#अन्यथा
      यदि (क्रमः == CblasColMajor) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, बफर, nthreads);
      } अन्यथा {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, बफर, nthreads);
      }
#यदि_अन्तः

  }
#यदि_अन्तः

  राशिः_FREE(बफर);

  FUNCTION_PROFILE_अन्तः(4, m * n + m + n, 2 * m * n);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
