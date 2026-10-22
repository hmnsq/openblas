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
#అన్తర్భావ <ctype.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XHBMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZHBMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CHBMV "
#అయితే_ముగింపు

స్థిర  పూర్ణాంకం (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, శూన్య *) = {
#నిర్వచితమైతే XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#లేదావేళ defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#లేకపోతే
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMPBUG
స్థిర  పూర్ణాంకం (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#లేదావేళ defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#లేకపోతే
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  అక్షరం uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  పూర్ణాంకం uplo;
  FLOAT *బఫరు;
#నిర్వచితమైతే SMPBUG
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;
  ఒకవేళ (uplo_arg  == 'V') uplo  = 2;
  ఒకవేళ (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  ఒకవేళ (incy == 0)          info = 11;
  ఒకవేళ (incx == 0)          info =  8;
  ఒకవేళ (lda  < k + 1)       info =  6;
  ఒకవేళ (k < 0)              info =  3;
  ఒకవేళ (n < 0)              info =  2;
  ఒకవేళ (uplo  < 0)          info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం,
	   గణన CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   శూన్య *VALPHA,
	   శూన్య  *va, blasint lda,
	   శూన్య  *vx, blasint incx,
	   శూన్య *VBETA,
	   శూన్య  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  FLOAT *బఫరు;
  పూర్ణాంకం uplo;
  blasint info;
#నిర్వచితమైతే SMPBUG
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (Uplo == Cblas_పై)         uplo  = 0;
    ఒకవేళ (Uplo == Cblas_క్రింది)         uplo  = 1;

    info = -1;

    ఒకవేళ (incy == 0)          info = 11;
    ఒకవేళ (incx == 0)          info =  8;
    ఒకవేళ (lda  < k + 1)       info =  6;
    ఒకవేళ (k < 0)              info =  3;
    ఒకవేళ (n < 0)              info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    ఒకవేళ (Uplo == Cblas_పై)         uplo  = 3;
    ఒకవేళ (Uplo == Cblas_క్రింది)         uplo  = 2;

    info = -1;

    ఒకవేళ (incy == 0)          info = 11;
    ఒకవేళ (incx == 0)          info =  8;
    ఒకవేళ (lda  < k + 1)       info =  6;
    ఒకవేళ (k < 0)              info =  3;
    ఒకవేళ (n < 0)              info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  ఒకవేళ (n == 0) తిరిగి;

  ఒకవేళ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ఖాళీ, 0, ఖాళీ, 0);

  ఒకవేళ ((alpha_r == ZERO) && (alpha_i == ZERO)) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  ఒకవేళ (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#నిర్వచితమైతే SMPBUG
  nthreads = అంకె_cpu_avail(2);

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, బఫరు);

#నిర్వచితమైతే SMPBUG
  } లేకపోతే {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(4, n * k / 2 + n,  n * k);

  IDEBUG_చివర;

  తిరిగి;
}
