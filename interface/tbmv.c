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
#define ERROR_NAME "QTBMV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DTBMV "
#अञ्ञथा
#define ERROR_NAME "STBMV "
#endif

स्थिर पूर्णाङ्क (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#नोचेद्यदि defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#अञ्ञथा
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#नोचेद्यदि defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#अञ्ञथा
  stbmv_thread_NUU, stbmv_thread_NUN, stbmv_thread_NLU, stbmv_thread_NLN,
  stbmv_thread_TUU, stbmv_thread_TUN, stbmv_thread_TLU, stbmv_thread_TLN,
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
  सचे (trans_arg == 'R') trans = 0;
  सचे (trans_arg == 'C') trans = 1;

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
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

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
    सचे (TransA == CblasConjNoTrans) trans = 0;
    सचे (TransA == CblasConjTrans)   trans = 1;

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
    सचे (TransA == CblasConjNoTrans) trans = 1;
    सचे (TransA == CblasConjTrans)   trans = 0;

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

  सचे (incx < 0 ) x -= (n - 1) * incx;

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

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  पटिनिवत्त;
}
