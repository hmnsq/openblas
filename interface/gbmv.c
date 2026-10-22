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

#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGBMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGBMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGBMV "
#అయితే_ముగింపు

స్థిర శూన్య (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, శూన్య *) = {
#నిర్వచితమైతే XDOUBLE
      qgbmv_n, qgbmv_t,
#లేదావేళ defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#లేకపోతే
      sgbmv_n, sgbmv_t,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#లేదావేళ defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#లేకపోతే
      sgbmv_thread_n, sgbmv_thread_t,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT
#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  అక్షరం trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *బఫరు;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ఒకవేళ (trans == 'N') i = 0;
  ఒకవేళ (trans == 'T') i = 1;
  ఒకవేళ (trans == 'R') i = 0;
  ఒకవేళ (trans == 'C') i = 1;

  ఒకవేళ (incy == 0)	 info = 13;
  ఒకవేళ (incx == 0)	 info = 10;
  ఒకవేళ (lda < kl + ku + 1) info = 8;
  ఒకవేళ (kl < 0)		 info = 5;
  ఒకవేళ (ku < 0)		 info = 4;
  ఒకవేళ (n < 0)		 info = 3;
  ఒకవేళ (m < 0)		 info = 2;
  ఒకవేళ (i < 0)		 info = 1;

  trans = i;

  ఒకవేళ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం,
	   గణన CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *బఫరు;
  blasint lenx, leny, info, t;
  పూర్ణాంకం trans;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (TransA == CblasNoTrans)     trans = 0;
    ఒకవేళ (TransA == CblasTrans)       trans = 1;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 0;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ఒకవేళ (incy == 0)	 info = 13;
    ఒకవేళ (incx == 0)	 info = 10;
    ఒకవేళ (lda < kl + ku + 1) info = 8;
    ఒకవేళ (kl < 0)		 info = 5;
    ఒకవేళ (ku < 0)		 info = 4;
    ఒకవేళ (n < 0)		 info = 3;
    ఒకవేళ (m < 0)		 info = 2;
    ఒకవేళ (trans < 0)	 info = 1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    ఒకవేళ (TransA == CblasNoTrans)     trans = 1;
    ఒకవేళ (TransA == CblasTrans)       trans = 0;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 1;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ఒకవేళ (incy == 0)	 info = 13;
    ఒకవేళ (incx == 0)	 info = 10;
    ఒకవేళ (lda < kl + ku + 1) info = 8;
    ఒకవేళ (kl < 0)		 info = 5;
    ఒకవేళ (ku < 0)		 info = 4;
    ఒకవేళ (n < 0)		 info = 3;
    ఒకవేళ (m < 0)		 info = 2;
    ఒకవేళ (trans < 0)	 info = 1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  ఒకవేళ ((m==0) || (n==0)) తిరిగి;

  lenx = n;
  leny = m;
  ఒకవేళ (trans) lenx = m;
  ఒకవేళ (trans) leny = n;

  ఒకవేళ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ఖాళీ, 0, ఖాళీ, 0);

  ఒకవేళ (alpha == ZERO) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0) x -= (lenx-1)*incx;
  ఒకవేళ (incy < 0) y -= (leny-1)*incy;

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#నిర్వచితమైతే SMP
  ఒకవేళ (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  లేకపోతే
    nthreads = అంకె_cpu_avail(2);

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  (gbmv[(పూర్ణాంకం)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {

    (gbmv_thread[(పూర్ణాంకం)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(1, m * n / 2 + n, m * n);

  IDEBUG_చివర;

  తిరిగి;
}
