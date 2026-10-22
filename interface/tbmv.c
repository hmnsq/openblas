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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTBMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DTBMV "
#ᬮᬶᬬ
#define ERROR_NAME "STBMV "
#endif

static puur_naa_nka (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#elif defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#ᬮᬶᬬ
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#endif
};

#ifdef SMP
static puur_naa_nka (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#elif defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#ᬮᬶᬬ
  stbmv_thread_NUU, stbmv_thread_NUN, stbmv_thread_NLU, stbmv_thread_NLN,
  stbmv_thread_TUU, stbmv_thread_TUN, stbmv_thread_TLU, stbmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANS, var_na *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;
  var_na diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  puur_naa_nka uplo;
  puur_naa_nka unit;
  puur_naa_nka trans;
  FLOAT *buffer;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ᬬᬾᬦ᭄ (trans_arg == 'N') trans = 0;
  ᬬᬾᬦ᭄ (trans_arg == 'T') trans = 1;
  ᬬᬾᬦ᭄ (trans_arg == 'R') trans = 0;
  ᬬᬾᬦ᭄ (trans_arg == 'C') trans = 1;

  ᬬᬾᬦ᭄ (diag_arg  == 'U') unit  = 0;
  ᬬᬾᬦ᭄ (diag_arg  == 'N') unit  = 1;

  ᬬᬾᬦ᭄ (uplo_arg  == 'U') uplo  = 0;
  ᬬᬾᬦ᭄ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᬬᬾᬦ᭄ (incx == 0)          info =  9;
  ᬬᬾᬦ᭄ (lda < k + 1)        info =  7;
  ᬬᬾᬦ᭄ (k < 0)              info =  5;
  ᬬᬾᬦ᭄ (n < 0)              info =  4;
  ᬬᬾᬦ᭄ (unit  < 0)          info =  3;
  ᬬᬾᬦ᭄ (trans < 0)          info =  2;
  ᬬᬾᬦ᭄ (uplo  < 0)          info =  1;

  ᬬᬾᬦ᭄ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᬯᬮᬶ;
  }

#ᬮᬶᬬ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  puur_naa_nka trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ᬬᬾᬦ᭄ (order == CblasColMajor) {
    ᬬᬾᬦ᭄ (Uplo == CblasUpper)         uplo  = 0;
    ᬬᬾᬦ᭄ (Uplo == CblasLower)         uplo  = 1;

    ᬬᬾᬦ᭄ (TransA == CblasNoTrans)     trans = 0;
    ᬬᬾᬦ᭄ (TransA == CblasTrans)       trans = 1;
    ᬬᬾᬦ᭄ (TransA == CblasConjNoTrans) trans = 0;
    ᬬᬾᬦ᭄ (TransA == CblasConjTrans)   trans = 1;

    ᬬᬾᬦ᭄ (Diag == CblasUnit)          unit  = 0;
    ᬬᬾᬦ᭄ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ᬬᬾᬦ᭄ (incx == 0)          info =  9;
    ᬬᬾᬦ᭄ (lda < k + 1)        info =  7;
    ᬬᬾᬦ᭄ (k < 0)              info =  5;
    ᬬᬾᬦ᭄ (n < 0)              info =  4;
    ᬬᬾᬦ᭄ (unit  < 0)          info =  3;
    ᬬᬾᬦ᭄ (trans < 0)          info =  2;
    ᬬᬾᬦ᭄ (uplo  < 0)          info =  1;
  }

  ᬬᬾᬦ᭄ (order == CblasRowMajor) {
    ᬬᬾᬦ᭄ (Uplo == CblasUpper)         uplo  = 1;
    ᬬᬾᬦ᭄ (Uplo == CblasLower)         uplo  = 0;

    ᬬᬾᬦ᭄ (TransA == CblasNoTrans)     trans = 1;
    ᬬᬾᬦ᭄ (TransA == CblasTrans)       trans = 0;
    ᬬᬾᬦ᭄ (TransA == CblasConjNoTrans) trans = 1;
    ᬬᬾᬦ᭄ (TransA == CblasConjTrans)   trans = 0;

    ᬬᬾᬦ᭄ (Diag == CblasUnit)          unit  = 0;
    ᬬᬾᬦ᭄ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ᬬᬾᬦ᭄ (incx == 0)          info =  9;
    ᬬᬾᬦ᭄ (lda < k + 1)        info =  7;
    ᬬᬾᬦ᭄ (k < 0)              info =  5;
    ᬬᬾᬦ᭄ (n < 0)              info =  4;
    ᬬᬾᬦ᭄ (unit  < 0)          info =  3;
    ᬬᬾᬦ᭄ (trans < 0)          info =  2;
    ᬬᬾᬦ᭄ (uplo  < 0)          info =  1;
  }

  ᬬᬾᬦ᭄ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᬯᬮᬶ;
  }

#endif

  ᬬᬾᬦ᭄ (n == 0) ᬯᬮᬶ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᬬᬾᬦ᭄ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᬬᬾᬦ᭄ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᬯᬮᬶ;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ᬬᬾᬦ᭄ (n < 50) nthreads = 1;
  ᬬᬾᬦ᭄ (nthreads > 2 && n < 500) nthreads = 2;

  ᬬᬾᬦ᭄ (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } ᬮᬶᬬ {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  ᬯᬮᬶ;
}
