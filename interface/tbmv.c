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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <ctype.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTBMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DTBMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "STBMV "
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#ಅನ್ಯಥಾ
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#endif
};

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#ಅನ್ಯಥಾ
  stbmv_thread_NUU, stbmv_thread_NUN, stbmv_thread_NLU, stbmv_thread_NLN,
  stbmv_thread_TUU, stbmv_thread_TUN, stbmv_thread_TLU, stbmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, ವರ್ಣ *TRANS, ವರ್ಣ *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  ವರ್ಣ uplo_arg  = *UPLO;
  ವರ್ಣ trans_arg = *TRANS;
  ವರ್ಣ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  ಪೂರ್ಣಾಙ್ಕ unit;
  ಪೂರ್ಣಾಙ್ಕ trans;
  FLOAT *buffer;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
  ಯದಿ (trans_arg == 'R') trans = 0;
  ಯದಿ (trans_arg == 'C') trans = 1;

  ಯದಿ (diag_arg  == 'U') unit  = 0;
  ಯದಿ (diag_arg  == 'N') unit  = 1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ಯದಿ (incx == 0)          info =  9;
  ಯದಿ (lda < k + 1)        info =  7;
  ಯದಿ (k < 0)              info =  5;
  ಯದಿ (n < 0)              info =  4;
  ಯದಿ (unit  < 0)          info =  3;
  ಯದಿ (trans < 0)          info =  2;
  ಯದಿ (uplo  < 0)          info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_UPLO Uplo,
	   ಗಣನಾ CBLAS_TRANSPOSE TransA, ಗಣನಾ CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  ಪೂರ್ಣಾಙ್ಕ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 0;
    ಯದಿ (Uplo == CblasLower)         uplo  = 1;

    ಯದಿ (TransA == CblasNoTrans)     trans = 0;
    ಯದಿ (TransA == CblasTrans)       trans = 1;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 0;
    ಯದಿ (TransA == CblasConjTrans)   trans = 1;

    ಯದಿ (Diag == CblasUnit)          unit  = 0;
    ಯದಿ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ಯದಿ (incx == 0)          info =  9;
    ಯದಿ (lda < k + 1)        info =  7;
    ಯದಿ (k < 0)              info =  5;
    ಯದಿ (n < 0)              info =  4;
    ಯದಿ (unit  < 0)          info =  3;
    ಯದಿ (trans < 0)          info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 1;
    ಯದಿ (Uplo == CblasLower)         uplo  = 0;

    ಯದಿ (TransA == CblasNoTrans)     trans = 1;
    ಯದಿ (TransA == CblasTrans)       trans = 0;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 1;
    ಯದಿ (TransA == CblasConjTrans)   trans = 0;

    ಯದಿ (Diag == CblasUnit)          unit  = 0;
    ಯದಿ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ಯದಿ (incx == 0)          info =  9;
    ಯದಿ (lda < k + 1)        info =  7;
    ಯದಿ (k < 0)              info =  5;
    ಯದಿ (n < 0)              info =  4;
    ಯದಿ (unit  < 0)          info =  3;
    ಯದಿ (trans < 0)          info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  ಯದಿ (n == 0) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ifdef SMP
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (n < 50) nthreads = 1;
  ಯದಿ (nthreads > 2 && n < 500) nthreads = 2;

  ಯದಿ (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } ಅನ್ಯಥಾ {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
