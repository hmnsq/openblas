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
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XTBMV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZTBMV "
#अञ्ञथा
#define ERROR_NAME "CTBMV "
#endif

स्थिर पूर्णाङ्क (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  xtbmv_NUU, xtbmv_NUN, xtbmv_NLU, xtbmv_NLN,
  xtbmv_TUU, xtbmv_TUN, xtbmv_TLU, xtbmv_TLN,
  xtbmv_RUU, xtbmv_RUN, xtbmv_RLU, xtbmv_RLN,
  xtbmv_CUU, xtbmv_CUN, xtbmv_CLU, xtbmv_CLN,
#नोचेद्यदि defined(DOUBLE)
  ztbmv_NUU, ztbmv_NUN, ztbmv_NLU, ztbmv_NLN,
  ztbmv_TUU, ztbmv_TUN, ztbmv_TLU, ztbmv_TLN,
  ztbmv_RUU, ztbmv_RUN, ztbmv_RLU, ztbmv_RLN,
  ztbmv_CUU, ztbmv_CUN, ztbmv_CLU, ztbmv_CLN,
#अञ्ञथा
  ctbmv_NUU, ctbmv_NUN, ctbmv_NLU, ctbmv_NLN,
  ctbmv_TUU, ctbmv_TUN, ctbmv_TLU, ctbmv_TLN,
  ctbmv_RUU, ctbmv_RUN, ctbmv_RLU, ctbmv_RLN,
  ctbmv_CUU, ctbmv_CUN, ctbmv_CLU, ctbmv_CLN,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xtbmv_thread_NUU, xtbmv_thread_NUN, xtbmv_thread_NLU, xtbmv_thread_NLN,
  xtbmv_thread_TUU, xtbmv_thread_TUN, xtbmv_thread_TLU, xtbmv_thread_TLN,
  xtbmv_thread_RUU, xtbmv_thread_RUN, xtbmv_thread_RLU, xtbmv_thread_RLN,
  xtbmv_thread_CUU, xtbmv_thread_CUN, xtbmv_thread_CLU, xtbmv_thread_CLN,
#नोचेद्यदि defined(DOUBLE)
  ztbmv_thread_NUU, ztbmv_thread_NUN, ztbmv_thread_NLU, ztbmv_thread_NLN,
  ztbmv_thread_TUU, ztbmv_thread_TUN, ztbmv_thread_TLU, ztbmv_thread_TLN,
  ztbmv_thread_RUU, ztbmv_thread_RUN, ztbmv_thread_RLU, ztbmv_thread_RLN,
  ztbmv_thread_CUU, ztbmv_thread_CUN, ztbmv_thread_CLU, ztbmv_thread_CLN,
#अञ्ञथा
  ctbmv_thread_NUU, ctbmv_thread_NUN, ctbmv_thread_NLU, ctbmv_thread_NLN,
  ctbmv_thread_TUU, ctbmv_thread_TUN, ctbmv_thread_TLU, ctbmv_thread_TLN,
  ctbmv_thread_RUU, ctbmv_thread_RUN, ctbmv_thread_RLU, ctbmv_thread_RLN,
  ctbmv_thread_CUU, ctbmv_thread_CUN, ctbmv_thread_CLU, ctbmv_thread_CLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS, वर्ण *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;
  वर्ण diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क unit;
  पूर्णाङ्क trans;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  सचे (trans_arg == 'N') trans = 0;
  सचे (trans_arg == 'T') trans = 1;
  सचे (trans_arg == 'R') trans = 2;
  सचे (trans_arg == 'C') trans = 3;

  सचे (diag_arg  == 'U') unit  = 0;
  सचे (diag_arg  == 'N') unit  = 1;

  सचे (uplo_arg  == 'U') uplo  = 0;
  सचे (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  सचे (incx == 0)          info =  9;
  सचे (lda < k + 1)        info =  7;
  सचे (k < 0)              info =  5;
  सचे (n < 0)              info =  4;
  सचे (unit  < 0)          info =  3;
  सचे (trans < 0)          info =  2;
  सचे (uplo  < 0)          info =  1;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_DIAG Diag,
	   blasint n, blasint k, शून्य  *va, blasint lda, शून्य  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;

  पूर्णाङ्क trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  सचे (order == CblasColMajor) {
    सचे (Uplo == CblasUpper)         uplo  = 0;
    सचे (Uplo == CblasLower)         uplo  = 1;

    सचे (TransA == CblasNoTrans)     trans = 0;
    सचे (TransA == CblasTrans)       trans = 1;
    सचे (TransA == CblasConjNoTrans) trans = 2;
    सचे (TransA == CblasConjTrans)   trans = 3;

    सचे (Diag == CblasUnit)          unit  = 0;
    सचे (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    सचे (incx == 0)          info =  9;
    सचे (lda < k + 1)        info =  7;
    सचे (k < 0)              info =  5;
    सचे (n < 0)              info =  4;
    सचे (unit  < 0)          info =  3;
    सचे (trans < 0)          info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (order == CblasRowMajor) {
    सचे (Uplo == CblasUpper)         uplo  = 1;
    सचे (Uplo == CblasLower)         uplo  = 0;

    सचे (TransA == CblasNoTrans)     trans = 1;
    सचे (TransA == CblasTrans)       trans = 0;
    सचे (TransA == CblasConjNoTrans) trans = 3;
    सचे (TransA == CblasConjTrans)   trans = 2;

    सचे (Diag == CblasUnit)          unit  = 0;
    सचे (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    सचे (incx == 0)          info =  9;
    सचे (lda < k + 1)        info =  7;
    सचे (k < 0)              info =  5;
    सचे (n < 0)              info =  4;
    सचे (unit  < 0)          info =  3;
    सचे (trans < 0)          info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  सचे (n == 0) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  सचे (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#ifdef SMP
  nthreads = अङ्कः_cpu_avail(2);

  सचे (n < 50) nthreads = 1;
  सचे (nthreads > 2 && n < 500) nthreads = 2;

  सचे (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } अञ्ञथा {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  पटिनिवत्त;
}
