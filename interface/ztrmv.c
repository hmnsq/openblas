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
#నిర్వచించు ERROR_NAME "XTRMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZTRMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CTRMV "
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#నిర్వచితమైతే XDOUBLE
  xtrmv_NUU, xtrmv_NUN, xtrmv_NLU, xtrmv_NLN,
  xtrmv_TUU, xtrmv_TUN, xtrmv_TLU, xtrmv_TLN,
  xtrmv_RUU, xtrmv_నిరంతరం, xtrmv_RLU, xtrmv_RLN,
  xtrmv_CUU, xtrmv_CUN, xtrmv_CLU, xtrmv_CLN,
#లేదావేళ defined(DOUBLE)
  ztrmv_NUU, ztrmv_NUN, ztrmv_NLU, ztrmv_NLN,
  ztrmv_TUU, ztrmv_TUN, ztrmv_TLU, ztrmv_TLN,
  ztrmv_RUU, ztrmv_నిరంతరం, ztrmv_RLU, ztrmv_RLN,
  ztrmv_CUU, ztrmv_CUN, ztrmv_CLU, ztrmv_CLN,
#లేకపోతే
  ctrmv_NUU, ctrmv_NUN, ctrmv_NLU, ctrmv_NLN,
  ctrmv_TUU, ctrmv_TUN, ctrmv_TLU, ctrmv_TLN,
  ctrmv_RUU, ctrmv_నిరంతరం, ctrmv_RLU, ctrmv_RLN,
  ctrmv_CUU, ctrmv_CUN, ctrmv_CLU, ctrmv_CLN,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  xtrmv_thread_NUU, xtrmv_thread_NUN, xtrmv_thread_NLU, xtrmv_thread_NLN,
  xtrmv_thread_TUU, xtrmv_thread_TUN, xtrmv_thread_TLU, xtrmv_thread_TLN,
  xtrmv_thread_RUU, xtrmv_thread_నిరంతరం, xtrmv_thread_RLU, xtrmv_thread_RLN,
  xtrmv_thread_CUU, xtrmv_thread_CUN, xtrmv_thread_CLU, xtrmv_thread_CLN,
#లేదావేళ defined(DOUBLE)
  ztrmv_thread_NUU, ztrmv_thread_NUN, ztrmv_thread_NLU, ztrmv_thread_NLN,
  ztrmv_thread_TUU, ztrmv_thread_TUN, ztrmv_thread_TLU, ztrmv_thread_TLN,
  ztrmv_thread_RUU, ztrmv_thread_నిరంతరం, ztrmv_thread_RLU, ztrmv_thread_RLN,
  ztrmv_thread_CUU, ztrmv_thread_CUN, ztrmv_thread_CLU, ztrmv_thread_CLN,
#లేకపోతే
  ctrmv_thread_NUU, ctrmv_thread_NUN, ctrmv_thread_NLU, ctrmv_thread_NLN,
  ctrmv_thread_TUU, ctrmv_thread_TUN, ctrmv_thread_TLU, ctrmv_thread_TLN,
  ctrmv_thread_RUU, ctrmv_thread_నిరంతరం, ctrmv_thread_RLU, ctrmv_thread_RLN,
  ctrmv_thread_CUU, ctrmv_thread_CUN, ctrmv_thread_CLU, ctrmv_thread_CLN,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, అక్షరం *TRANS, అక్షరం *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  అక్షరం uplo_arg  = *UPLO;
  అక్షరం trans_arg = *TRANS;
  అక్షరం diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  పూర్ణాంకం uplo;
  పూర్ణాంకం unit;
  పూర్ణాంకం trans, బఫరు_కొలత;
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

  ఒకవేళ (incx == 0)          info =  8;
  ఒకవేళ (lda  < MAX(1, n))   info =  6;
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
	   blasint n, శూన్య  *va, blasint lda, శూన్య  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  పూర్ణాంకం trans, uplo, unit, బఫరు_కొలత;
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
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 2;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 3;

    ఒకవేళ (Diag == CblasUnit)          unit  = 0;
    ఒకవేళ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ఒకవేళ (incx == 0)          info =  8;
    ఒకవేళ (lda  < MAX(1, n))   info =  6;
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

    ఒకవేళ (incx == 0)          info =  8;
    ఒకవేళ (lda  < MAX(1, n))   info =  6;
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

#నిర్వచితమైతే SMP
  // Calibrated on a Xeon E5-2630
  ఒకవేళ(1L * n * n > 36L * పరిమాణం(FLOAT) * పరిమాణం(FLOAT) * GEMM_MULTITHREAD_పరిమితి) {
    nthreads = అంకె_cpu_avail(2);
    ఒకవేళ(nthreads > 2 && 1L * n * n < 64L * పరిమాణం(FLOAT) * పరిమాణం(FLOAT) * GEMM_MULTITHREAD_పరిమితి)
      nthreads = 2;
  } లేకపోతే
      nthreads = 1;

  ఒకవేళ(nthreads > 1) {
    బఫరు_కొలత = n > 16 ? 0 : n * 4 + 40;
  }
  లేకపోతే
#అయితే_ముగింపు
  {
    బఫరు_కొలత = ((n - 1) / DTB_ENTRIES) * 2 * DTB_ENTRIES + 32 / పరిమాణం(FLOAT);
    // It seems to be required for some K8 or Barcelona CPU
    బఫరు_కొలత += 8;
    ఒకవేళ(incx != 1)
      బఫరు_కొలత += n * 2;
  }
  దొంతర_ALLOC(బఫరు_కొలత, FLOAT, బఫరు);

#నిర్వచితమైతే SMP
  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  దొంతర_FREE(బఫరు);

  FUNCTION_PROFILE_చివర(4, n * n / 2 + n,  n * n);

  IDEBUG_చివర;

  తిరిగి;
}
