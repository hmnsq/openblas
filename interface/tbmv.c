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
#నిర్వచించు ERROR_NAME "QTBMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DTBMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "STBMV "
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, శూన్య *) = {
#నిర్వచితమైతే XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#లేదావేళ defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#లేకపోతే
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#లేదావేళ defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#లేకపోతే
  stbmv_thread_NUU, stbmv_thread_NUN, stbmv_thread_NLU, stbmv_thread_NLN,
  stbmv_thread_TUU, stbmv_thread_TUN, stbmv_thread_TLU, stbmv_thread_TLN,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, అక్షరం *TRANS, అక్షరం *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  అక్షరం uplo_arg  = *UPLO;
  అక్షరం trans_arg = *TRANS;
  అక్షరం diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  పూర్ణాంకం uplo;
  పూర్ణాంకం unit;
  పూర్ణాంకం trans;
  FLOAT *బఫరు;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'T') trans = 1;
  ఒకవేళ (trans_arg == 'R') trans = 0;
  ఒకవేళ (trans_arg == 'C') trans = 1;

  ఒకవేళ (diag_arg  == 'U') unit  = 0;
  ఒకవేళ (diag_arg  == 'N') unit  = 1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ఒకవేళ (incx == 0)          info =  9;
  ఒకవేళ (lda < k + 1)        info =  7;
  ఒకవేళ (k < 0)              info =  5;
  ఒకవేళ (n < 0)              info =  4;
  ఒకవేళ (unit  < 0)          info =  3;
  ఒకవేళ (trans < 0)          info =  2;
  ఒకవేళ (uplo  < 0)          info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_UPLO Uplo,
	   గణన CBLAS_TRANSPOSE TransA, గణన CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  పూర్ణాంకం trans, uplo, unit;
  blasint info;
  FLOAT *బఫరు;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (Uplo == Cblas_పై)         uplo  = 0;
    ఒకవేళ (Uplo == Cblas_క్రింది)         uplo  = 1;

    ఒకవేళ (TransA == CblasNoTrans)     trans = 0;
    ఒకవేళ (TransA == CblasTrans)       trans = 1;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 0;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 1;

    ఒకవేళ (Diag == CblasUnit)          unit  = 0;
    ఒకవేళ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ఒకవేళ (incx == 0)          info =  9;
    ఒకవేళ (lda < k + 1)        info =  7;
    ఒకవేళ (k < 0)              info =  5;
    ఒకవేళ (n < 0)              info =  4;
    ఒకవేళ (unit  < 0)          info =  3;
    ఒకవేళ (trans < 0)          info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    ఒకవేళ (Uplo == Cblas_పై)         uplo  = 1;
    ఒకవేళ (Uplo == Cblas_క్రింది)         uplo  = 0;

    ఒకవేళ (TransA == CblasNoTrans)     trans = 1;
    ఒకవేళ (TransA == CblasTrans)       trans = 0;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 1;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 0;

    ఒకవేళ (Diag == CblasUnit)          unit  = 0;
    ఒకవేళ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ఒకవేళ (incx == 0)          info =  9;
    ఒకవేళ (lda < k + 1)        info =  7;
    ఒకవేళ (k < 0)              info =  5;
    ఒకవేళ (n < 0)              info =  4;
    ఒకవేళ (unit  < 0)          info =  3;
    ఒకవేళ (trans < 0)          info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  ఒకవేళ (n == 0) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0 ) x -= (n - 1) * incx;

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#నిర్వచితమైతే SMP
  nthreads = అంకె_cpu_avail(2);

  ఒకవేళ (n < 50) nthreads = 1;
  ఒకవేళ (nthreads > 2 && n < 500) nthreads = 2;

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(1, n * k / 2 + n, n * k);

  IDEBUG_చివర;

  తిరిగి;
}
