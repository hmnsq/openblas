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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
#నిర్వచితమైతే __64BIT__
#నిర్వచించు SMPTEST 1
#అయితే_ముగింపు
#అయితే_ముగింపు


#నిర్వచితమైతే XDOUBLE
#నిర్వచితం_కాకపోతే CONJ
#నిర్వచించు ERROR_NAME "XGERU  "
#లేకపోతే
#నిర్వచించు ERROR_NAME "XGERC  "
#అయితే_ముగింపు
#లేదావేళ defined DOUBLE
#నిర్వచితం_కాకపోతే CONJ
#నిర్వచించు ERROR_NAME "ZGERU  "
#లేకపోతే
#నిర్వచించు ERROR_NAME "ZGERC  "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితం_కాకపోతే CONJ
#నిర్వచించు ERROR_NAME "CGERU  "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGERC "
#అయితే_ముగింపు
#అయితే_ముగింపు

#ఒకవేళ   defined XDOUBLE
#నిర్వచితం_కాకపోతే CONJ
#నిర్వచించు GER		GERU_K
#నిర్వచించు GER_THREAD	xger_thread_U
#లేకపోతే
#నిర్వచించు GER		GERC_K
#నిర్వచించు GER_THREAD	xger_thread_C
#నిర్వచించు GERV		GERV_K
#నిర్వచించు GERV_THREAD	xger_thread_V
#అయితే_ముగింపు
#లేదావేళ defined DOUBLE
#నిర్వచితం_కాకపోతే CONJ
#నిర్వచించు GER		GERU_K
#నిర్వచించు GER_THREAD	zger_thread_U
#లేకపోతే
#నిర్వచించు GER		GERC_K
#నిర్వచించు GER_THREAD	zger_thread_C
#నిర్వచించు GERV		GERV_K
#నిర్వచించు GERV_THREAD	zger_thread_V
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితం_కాకపోతే CONJ
#నిర్వచించు GER		GERU_K
#నిర్వచించు GER_THREAD	cger_thread_U
#లేకపోతే
#నిర్వచించు GER		GERC_K
#నిర్వచించు GER_THREAD	cger_thread_C
#నిర్వచించు GERV		GERV_K
#నిర్వచించు GERV_THREAD	cger_thread_V
#అయితే_ముగింపు
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(blasint *M, blasint *N, FLOAT *Alpha,
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
  FLOAT *బఫరు;
#నిర్వచితమైతే SMPTEST
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  ఒకవేళ (lda < MAX(1,m)) info = 9;
  ఒకవేళ (incy == 0)      info = 7;
  ఒకవేళ (incx == 0)      info = 5;
  ఒకవేళ (n < 0)          info = 2;
  ఒకవేళ (m < 0)          info = 1;

  ఒకవేళ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం,
	   blasint m, blasint n,
	   శూన్య *VAlpha,
	   శూన్య  *vx, blasint incx,
	   శూన్య  *vy, blasint incy,
	   శూన్య  *va, blasint lda) {

  FLOAT* Alpha = (FLOAT*) VAlpha;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT  alpha_r = Alpha[0];
  FLOAT  alpha_i = Alpha[1];

  FLOAT *బఫరు;
  blasint info, t;
#నిర్వచితమైతే SMPTEST
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    info = -1;

    ఒకవేళ (lda < MAX(1,m)) info = 9;
    ఒకవేళ (incy == 0)      info = 7;
    ఒకవేళ (incx == 0)      info = 5;
    ఒకవేళ (n < 0)          info = 2;
    ఒకవేళ (m < 0)          info = 1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    బఫరు = x;
    x = y;
    y = బఫరు;

    ఒకవేళ (lda < MAX(1,m)) info = 9;
    ఒకవేళ (incy == 0)      info = 7;
    ఒకవేళ (incx == 0)      info = 5;
    ఒకవేళ (n < 0)          info = 2;
    ఒకవేళ (m < 0)          info = 1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  /*     Quick return if possible. */
  ఒకవేళ (m == 0 || n == 0) తిరిగి;

  ఒకవేళ ((alpha_r == 0.) && (alpha_i == 0.)) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incy < 0) y -= (n - 1) * incy * 2;
  ఒకవేళ (incx < 0) x -= (m - 1) * incx * 2;

  దొంతర_ALLOC(2 * m, FLOAT, బఫరు);

#నిర్వచితమైతే SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  ఒకవేళ(1L * m * n > 36L * పరిమాణం(FLOAT) * పరిమాణం(FLOAT) * GEMM_MULTITHREAD_పరిమితి)
    nthreads = అంకె_cpu_avail(2);
  లేకపోతే
    nthreads = 1;

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

#ఒకవేళ !defined(CBLAS) || !defined(CONJ)
  GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, బఫరు);
#లేకపోతే
  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    GER(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, బఫరు);
  } లేకపోతే {
    GERV(m, n, 0, alpha_r, alpha_i, x, incx, y, incy, a, lda, బఫరు);
  }
#అయితే_ముగింపు

#నిర్వచితమైతే SMPTEST

  } లేకపోతే {

#ఒకవేళ !defined(CBLAS) || !defined(CONJ)
      GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, బఫరు, nthreads);
#లేకపోతే
      ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
	GER_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, బఫరు, nthreads);
      } లేకపోతే {
	GERV_THREAD(m, n, Alpha, x, incx, y, incy, a, lda, బఫరు, nthreads);
      }
#అయితే_ముగింపు

  }
#అయితే_ముగింపు

  దొంతర_FREE(బఫరు);

  FUNCTION_PROFILE_చివర(4, m * n + m + n, 2 * m * n);

  IDEBUG_చివర;

  తిరిగి;

}
