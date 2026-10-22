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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <ctype.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XTBMV "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZTBMV "
#ଅନ୍ୟଥା
#define ERROR_NAME "CTBMV "
#endif

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ଶୂନ୍ୟ *) = {
#ifdef XDOUBLE
  xtbmv_NUU, xtbmv_NUN, xtbmv_NLU, xtbmv_NLN,
  xtbmv_TUU, xtbmv_TUN, xtbmv_TLU, xtbmv_TLN,
  xtbmv_RUU, xtbmv_RUN, xtbmv_RLU, xtbmv_RLN,
  xtbmv_CUU, xtbmv_CUN, xtbmv_CLU, xtbmv_CLN,
#ନଚେତଯଦି defined(DOUBLE)
  ztbmv_NUU, ztbmv_NUN, ztbmv_NLU, ztbmv_NLN,
  ztbmv_TUU, ztbmv_TUN, ztbmv_TLU, ztbmv_TLN,
  ztbmv_RUU, ztbmv_RUN, ztbmv_RLU, ztbmv_RLN,
  ztbmv_CUU, ztbmv_CUN, ztbmv_CLU, ztbmv_CLN,
#ଅନ୍ୟଥା
  ctbmv_NUU, ctbmv_NUN, ctbmv_NLU, ctbmv_NLN,
  ctbmv_TUU, ctbmv_TUN, ctbmv_TLU, ctbmv_TLN,
  ctbmv_RUU, ctbmv_RUN, ctbmv_RLU, ctbmv_RLN,
  ctbmv_CUU, ctbmv_CUN, ctbmv_CLU, ctbmv_CLN,
#endif
};

#ifdef SMP
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ପୂର୍ଣ୍ଣସଂଖ୍ୟା) = {
#ifdef XDOUBLE
  xtbmv_thread_NUU, xtbmv_thread_NUN, xtbmv_thread_NLU, xtbmv_thread_NLN,
  xtbmv_thread_TUU, xtbmv_thread_TUN, xtbmv_thread_TLU, xtbmv_thread_TLN,
  xtbmv_thread_RUU, xtbmv_thread_RUN, xtbmv_thread_RLU, xtbmv_thread_RLN,
  xtbmv_thread_CUU, xtbmv_thread_CUN, xtbmv_thread_CLU, xtbmv_thread_CLN,
#ନଚେତଯଦି defined(DOUBLE)
  ztbmv_thread_NUU, ztbmv_thread_NUN, ztbmv_thread_NLU, ztbmv_thread_NLN,
  ztbmv_thread_TUU, ztbmv_thread_TUN, ztbmv_thread_TLU, ztbmv_thread_TLN,
  ztbmv_thread_RUU, ztbmv_thread_RUN, ztbmv_thread_RLU, ztbmv_thread_RLN,
  ztbmv_thread_CUU, ztbmv_thread_CUN, ztbmv_thread_CLU, ztbmv_thread_CLN,
#ଅନ୍ୟଥା
  ctbmv_thread_NUU, ctbmv_thread_NUN, ctbmv_thread_NLU, ctbmv_thread_NLN,
  ctbmv_thread_TUU, ctbmv_thread_TUN, ctbmv_thread_TLU, ctbmv_thread_TLN,
  ctbmv_thread_RUU, ctbmv_thread_RUN, ctbmv_thread_RLU, ctbmv_thread_RLN,
  ctbmv_thread_CUU, ctbmv_thread_CUN, ctbmv_thread_CLU, ctbmv_thread_CLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, ଅକ୍ଷର *TRANS, ଅକ୍ଷର *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  ଅକ୍ଷର uplo_arg  = *UPLO;
  ଅକ୍ଷର trans_arg = *TRANS;
  ଅକ୍ଷର diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା unit;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans;
  FLOAT *buffer;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'T') trans = 1;
  ଯଦି (trans_arg == 'R') trans = 2;
  ଯଦି (trans_arg == 'C') trans = 3;

  ଯଦି (diag_arg  == 'U') unit  = 0;
  ଯଦି (diag_arg  == 'N') unit  = 1;

  ଯଦି (uplo_arg  == 'U') uplo  = 0;
  ଯଦି (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ଯଦି (incx == 0)          info =  9;
  ଯଦି (lda < k + 1)        info =  7;
  ଯଦି (k < 0)              info =  5;
  ଯଦି (n < 0)              info =  4;
  ଯଦି (unit  < 0)          info =  3;
  ଯଦି (trans < 0)          info =  2;
  ଯଦି (uplo  < 0)          info =  1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_UPLO Uplo,
	   ଗଣନା CBLAS_TRANSPOSE TransA, ଗଣନା CBLAS_DIAG Diag,
	   blasint n, blasint k, ଶୂନ୍ୟ  *va, blasint lda, ଶୂନ୍ୟ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (Uplo == CblasUpper)         uplo  = 0;
    ଯଦି (Uplo == CblasLower)         uplo  = 1;

    ଯଦି (TransA == CblasNoTrans)     trans = 0;
    ଯଦି (TransA == CblasTrans)       trans = 1;
    ଯଦି (TransA == CblasConjNoTrans) trans = 2;
    ଯଦି (TransA == CblasConjTrans)   trans = 3;

    ଯଦି (Diag == CblasUnit)          unit  = 0;
    ଯଦି (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ଯଦି (incx == 0)          info =  9;
    ଯଦି (lda < k + 1)        info =  7;
    ଯଦି (k < 0)              info =  5;
    ଯଦି (n < 0)              info =  4;
    ଯଦି (unit  < 0)          info =  3;
    ଯଦି (trans < 0)          info =  2;
    ଯଦି (uplo  < 0)          info =  1;
  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (Uplo == CblasUpper)         uplo  = 1;
    ଯଦି (Uplo == CblasLower)         uplo  = 0;

    ଯଦି (TransA == CblasNoTrans)     trans = 1;
    ଯଦି (TransA == CblasTrans)       trans = 0;
    ଯଦି (TransA == CblasConjNoTrans) trans = 3;
    ଯଦି (TransA == CblasConjTrans)   trans = 2;

    ଯଦି (Diag == CblasUnit)          unit  = 0;
    ଯଦି (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ଯଦି (incx == 0)          info =  9;
    ଯଦି (lda < k + 1)        info =  7;
    ଯଦି (k < 0)              info =  5;
    ଯଦି (n < 0)              info =  4;
    ଯଦି (unit  < 0)          info =  3;
    ଯଦି (trans < 0)          info =  2;
    ଯଦି (uplo  < 0)          info =  1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#endif

  ଯଦି (n == 0) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ଯଦି (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ifdef SMP
  nthreads = ଅଙ୍କ_cpu_avail(2);

  ଯଦି (n < 50) nthreads = 1;
  ଯଦି (nthreads > 2 && n < 500) nthreads = 2;

  ଯଦି (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } ଅନ୍ୟଥା {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
