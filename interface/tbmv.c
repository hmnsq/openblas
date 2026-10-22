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
#නැතහොත්
#define ERROR_NAME "STBMV "
#endif

static පූර්ණ_සංඛ්‍යා (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#elif defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#නැතහොත්
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#endif
};

#ifdef SMP
static පූර්ණ_සංඛ්‍යා (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, පූර්ණ_සංඛ්‍යා) = {
#ifdef XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#elif defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#නැතහොත්
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
  පූර්ණ_සංඛ්‍යා uplo;
  පූර්ණ_සංඛ්‍යා unit;
  පූර්ණ_සංඛ්‍යා trans;
  FLOAT *buffer;
#ifdef SMP
  පූර්ණ_සංඛ්‍යා nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  නම් (trans_arg == 'N') trans = 0;
  නම් (trans_arg == 'T') trans = 1;
  නම් (trans_arg == 'R') trans = 0;
  නම් (trans_arg == 'C') trans = 1;

  නම් (diag_arg  == 'U') unit  = 0;
  නම් (diag_arg  == 'N') unit  = 1;

  නම් (uplo_arg  == 'U') uplo  = 0;
  නම් (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  නම් (incx == 0)          info =  9;
  නම් (lda < k + 1)        info =  7;
  නම් (k < 0)              info =  5;
  නම් (n < 0)              info =  4;
  නම් (unit  < 0)          info =  3;
  නම් (trans < 0)          info =  2;
  නම් (uplo  < 0)          info =  1;

  නම් (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ආපසු;
  }

#නැතහොත්

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  පූර්ණ_සංඛ්‍යා trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  පූර්ණ_සංඛ්‍යා nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  නම් (order == CblasColMajor) {
    නම් (Uplo == CblasUpper)         uplo  = 0;
    නම් (Uplo == CblasLower)         uplo  = 1;

    නම් (TransA == CblasNoTrans)     trans = 0;
    නම් (TransA == CblasTrans)       trans = 1;
    නම් (TransA == CblasConjNoTrans) trans = 0;
    නම් (TransA == CblasConjTrans)   trans = 1;

    නම් (Diag == CblasUnit)          unit  = 0;
    නම් (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    නම් (incx == 0)          info =  9;
    නම් (lda < k + 1)        info =  7;
    නම් (k < 0)              info =  5;
    නම් (n < 0)              info =  4;
    නම් (unit  < 0)          info =  3;
    නම් (trans < 0)          info =  2;
    නම් (uplo  < 0)          info =  1;
  }

  නම් (order == CblasRowMajor) {
    නම් (Uplo == CblasUpper)         uplo  = 1;
    නම් (Uplo == CblasLower)         uplo  = 0;

    නම් (TransA == CblasNoTrans)     trans = 1;
    නම් (TransA == CblasTrans)       trans = 0;
    නම් (TransA == CblasConjNoTrans) trans = 1;
    නම් (TransA == CblasConjTrans)   trans = 0;

    නම් (Diag == CblasUnit)          unit  = 0;
    නම් (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    නම් (incx == 0)          info =  9;
    නම් (lda < k + 1)        info =  7;
    නම් (k < 0)              info =  5;
    නම් (n < 0)              info =  4;
    නම් (unit  < 0)          info =  3;
    නම් (trans < 0)          info =  2;
    නම් (uplo  < 0)          info =  1;
  }

  නම් (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ආපසු;
  }

#endif

  නම් (n == 0) ආපසු;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  නම් (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  නම් (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ආපසු;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  නම් (n < 50) nthreads = 1;
  නම් (nthreads > 2 && n < 500) nthreads = 2;

  නම් (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } නැතහොත් {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  ආපසු;
}
