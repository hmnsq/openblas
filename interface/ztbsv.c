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
#నిర్వచించు ERROR_NAME "XTBSV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZTBSV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CTBSV "
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, శూన్య *) = {
#నిర్వచితమైతే XDOUBLE
  xtbsv_NUU, xtbsv_NUN, xtbsv_NLU, xtbsv_NLN,
  xtbsv_TUU, xtbsv_TUN, xtbsv_TLU, xtbsv_TLN,
  xtbsv_RUU, xtbsv_నిరంతరం, xtbsv_RLU, xtbsv_RLN,
  xtbsv_CUU, xtbsv_CUN, xtbsv_CLU, xtbsv_CLN,
#లేదావేళ defined(DOUBLE)
  ztbsv_NUU, ztbsv_NUN, ztbsv_NLU, ztbsv_NLN,
  ztbsv_TUU, ztbsv_TUN, ztbsv_TLU, ztbsv_TLN,
  ztbsv_RUU, ztbsv_నిరంతరం, ztbsv_RLU, ztbsv_RLN,
  ztbsv_CUU, ztbsv_CUN, ztbsv_CLU, ztbsv_CLN,
#లేకపోతే
  ctbsv_NUU, ctbsv_NUN, ctbsv_NLU, ctbsv_NLN,
  ctbsv_TUU, ctbsv_TUN, ctbsv_TLU, ctbsv_TLN,
  ctbsv_RUU, ctbsv_నిరంతరం, ctbsv_RLU, ctbsv_RLN,
  ctbsv_CUU, ctbsv_CUN, ctbsv_CLU, ctbsv_CLN,
#అయితే_ముగింపు
};

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

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'T') trans = 1;
  ఒకవేళ (trans_arg == 'R') trans = 2;
  ఒకవేళ (trans_arg == 'C') trans = 3;

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
	   blasint n, blasint k, శూన్య  *va, blasint lda, శూన్య  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  పూర్ణాంకం trans, uplo, unit;
  blasint info;
  FLOAT *బఫరు;

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
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 2;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 3;

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
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 3;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 2;

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

  ఒకవేళ (incx < 0 ) x -= (n - 1) * incx * 2;

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, బఫరు);

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(4, n * k / 2 + n, n * k);

  IDEBUG_చివర;

  తిరిగి;
}
