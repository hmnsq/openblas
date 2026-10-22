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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित SMP
#अगर_घोषित __64BIT__
#घोषणा SMPTEST 1
#पूर्ण_अगर
#पूर्ण_अगर


#अगर_घोषित XDOUBLE
#अगर_अघोषित CONJ
#घोषणा ERROR_NAME "XGERU  "
#अन्यथा
#घोषणा ERROR_NAME "XGERC  "
#पूर्ण_अगर
#या_अगर defined DOUBLE
#अगर_अघोषित CONJ
#घोषणा ERROR_NAME "ZGERU  "
#अन्यथा
#घोषणा ERROR_NAME "ZGERC  "
#पूर्ण_अगर
#अन्यथा
#अगर_अघोषित CONJ
#घोषणा ERROR_NAME "CGERU  "
#अन्यथा
#घोषणा ERROR_NAME "CGERC "
#पूर्ण_अगर
#पूर्ण_अगर

#अगर   defined XDOUBLE
#अगर_अघोषित CONJ
#घोषणा GER		GERU_K
#घोषणा GER_THREAD	xger_thread_U
#अन्यथा
#घोषणा GER		GERC_K
#घोषणा GER_THREAD	xger_thread_C
#घोषणा GERV		GERV_K
#घोषणा GERV_THREAD	xger_thread_V
#पूर्ण_अगर
#या_अगर defined DOUBLE
#अगर_अघोषित CONJ
#घोषणा GER		GERU_K
#घोषणा GER_THREAD	zger_thread_U
#अन्यथा
#घोषणा GER		GERC_K
#घोषणा GER_THREAD	zger_thread_C
#घोषणा GERV		GERV_K
#घोषणा GERV_THREAD	zger_thread_V
#पूर्ण_अगर
#अन्यथा
#अगर_अघोषित CONJ
#घोषणा GER		GERU_K
#घोषणा GER_THREAD	cger_thread_U
#अन्यथा
#घोषणा GER		GERC_K
#घोषणा GER_THREAD	cger_thread_C
#घोषणा GERV		GERV_K
#घोषणा GERV_THREAD	cger_thread_V
#पूर्ण_अगर
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
#अगर_घोषित SMPTEST
  पूर्णांक nthreads;
#पूर्ण_अगर

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  अगर (lda < MAX(1,m)) info = 9;
  अगर (incy == 0)      info = 7;
  अगर (incx == 0)      info = 5;
  अगर (n < 0)          info = 2;
  अगर (m < 0)          info = 1;

  अगर (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   blasint m, blasint n,
	   व्योम *VAlpha,
	   व्योम  *vx, blasint incx,
	   व्योम  *vy, blasint incy,
	   व्योम  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *बफर;
  blasint info, t;
#अगर_घोषित SMPTEST
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  info  =  0;

  अगर (order == CblasColMajor) {
    info = -1;

    अगर (lda < MAX(1,m)) info = 9;
    अगर (incy == 0)      info = 7;
    अगर (incx == 0)      info = 5;
    अगर (n < 0)          info = 2;
    अगर (m < 0)          info = 1;
  }

  अगर (order == CblasRowMajor) {
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

    अगर (lda < MAX(1,m)) info = 9;
    अगर (incy == 0)      info = 7;
    अगर (incx == 0)      info = 5;
    अगर (n < 0)          info = 2;
    अगर (m < 0)          info = 1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर

  /*     Quick return if possible. */
  अगर (m == 0 || n == 0) वापस;

  अगर ((alpha_r == 0.) && (alpha_i == 0.)) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incy < 0) y -= (n - 1) * incy * 2;
  अगर (incx < 0) x -= (m - 1) * incx * 2;

  ढेर_ALLOC(2 * m, FLOAT, बफर);

#अगर_घोषित SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  अगर(1L * m * n > 36L * माप(FLOAT) * माप(FLOAT) * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = अंक_cpu_avail(2);
  अन्यथा
    nthreads = 1;

  अगर (nthreads == 1) {
#पूर्ण_अगर

#अगर !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, बफर);
#अन्यथा
  अगर (order == CblasColMajor) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, बफर);
  } अन्यथा {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, बफर);
  }
#पूर्ण_अगर

#अगर_घोषित SMPTEST

  } अन्यथा {

#अगर !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, बफर, nthreads);
#अन्यथा
      अगर (order == CblasColMajor) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, बफर, nthreads);
      } अन्यथा {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, बफर, nthreads);
      }
#पूर्ण_अगर

  }
#पूर्ण_अगर

  ढेर_FREE(बफर);

  FUNCTION_PROFILE_END(4, m * n + m + n, 2 * m * n);

  IDEBUG_END;

  वापस;

}
