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
#నిర్వచించు ERROR_NAME "XTPMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZTPMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CTPMV "
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, శూన్య *) = {
#నిర్వచితమైతే XDOUBLE
  xtpmv_NUU, xtpmv_NUN, xtpmv_NLU, xtpmv_NLN,
  xtpmv_TUU, xtpmv_TUN, xtpmv_TLU, xtpmv_TLN,
  xtpmv_RUU, xtpmv_నిరంతరం, xtpmv_RLU, xtpmv_RLN,
  xtpmv_CUU, xtpmv_CUN, xtpmv_CLU, xtpmv_CLN,
#లేదావేళ defined(DOUBLE)
  ztpmv_NUU, ztpmv_NUN, ztpmv_NLU, ztpmv_NLN,
  ztpmv_TUU, ztpmv_TUN, ztpmv_TLU, ztpmv_TLN,
  ztpmv_RUU, ztpmv_నిరంతరం, ztpmv_RLU, ztpmv_RLN,
  ztpmv_CUU, ztpmv_CUN, ztpmv_CLU, ztpmv_CLN,
#లేకపోతే
  ctpmv_NUU, ctpmv_NUN, ctpmv_NLU, ctpmv_NLN,
  ctpmv_TUU, ctpmv_TUN, ctpmv_TLU, ctpmv_TLN,
  ctpmv_RUU, ctpmv_నిరంతరం, ctpmv_RLU, ctpmv_RLN,
  ctpmv_CUU, ctpmv_CUN, ctpmv_CLU, ctpmv_CLN,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  xtpmv_thread_NUU, xtpmv_thread_NUN, xtpmv_thread_NLU, xtpmv_thread_NLN,
  xtpmv_thread_TUU, xtpmv_thread_TUN, xtpmv_thread_TLU, xtpmv_thread_TLN,
  xtpmv_thread_RUU, xtpmv_thread_నిరంతరం, xtpmv_thread_RLU, xtpmv_thread_RLN,
  xtpmv_thread_CUU, xtpmv_thread_CUN, xtpmv_thread_CLU, xtpmv_thread_CLN,
#లేదావేళ defined(DOUBLE)
  ztpmv_thread_NUU, ztpmv_thread_NUN, ztpmv_thread_NLU, ztpmv_thread_NLN,
  ztpmv_thread_TUU, ztpmv_thread_TUN, ztpmv_thread_TLU, ztpmv_thread_TLN,
  ztpmv_thread_RUU, ztpmv_thread_నిరంతరం, ztpmv_thread_RLU, ztpmv_thread_RLN,
  ztpmv_thread_CUU, ztpmv_thread_CUN, ztpmv_thread_CLU, ztpmv_thread_CLN,
#లేకపోతే
  ctpmv_thread_NUU, ctpmv_thread_NUN, ctpmv_thread_NLU, ctpmv_thread_NLN,
  ctpmv_thread_TUU, ctpmv_thread_TUN, ctpmv_thread_TLU, ctpmv_thread_TLN,
  ctpmv_thread_RUU, ctpmv_thread_నిరంతరం, ctpmv_thread_RLU, ctpmv_thread_RLN,
  ctpmv_thread_CUU, ctpmv_thread_CUN, ctpmv_thread_CLU, ctpmv_thread_CLN,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, అక్షరం *TRANS, అక్షరం *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  అక్షరం uplo_arg  = *UPLO;
  అక్షరం trans_arg = *TRANS;
  అక్షరం diag_arg  = *DIAG;

  blasint n    = *N;
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
  ఒకవేళ (trans_arg == 'R') trans = 2;
  ఒకవేళ (trans_arg == 'C') trans = 3;

  ఒకవేళ (diag_arg  == 'U') unit  = 0;
  ఒకవేళ (diag_arg  == 'N') unit  = 1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ఒకవేళ (incx == 0)          info =  7;
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
	   blasint n, శూన్య  *va, శూన్య  *vx, blasint incx) {

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
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

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

    ఒకవేళ (incx == 0)          info =  7;
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

    ఒకవేళ (incx == 0)          info =  7;
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

#నిర్వచితమైతే SMP
  nthreads = అంకె_cpu_avail(2);

  ఒకవేళ (n < 50) nthreads = 1;
  ఒకవేళ (nthreads > 2 && n < 500) nthreads = 2;

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, బఫరు);

#నిర్వచితమైతే SMP

  } లేకపోతే {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(4, n * n / 2 + n, n * n);

  IDEBUG_చివర;

  తిరిగి;
}
