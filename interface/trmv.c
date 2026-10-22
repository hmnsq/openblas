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
#define ERROR_NAME "QTRMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DTRMV "
#නැතහොත්
#define ERROR_NAME "STRMV "
#endif

static පූර්ණ_සංඛ්‍යා (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  qtrmv_NUU, qtrmv_NUN, qtrmv_NLU, qtrmv_NLN,
  qtrmv_TUU, qtrmv_TUN, qtrmv_TLU, qtrmv_TLN,
#elif defined(DOUBLE)
  dtrmv_NUU, dtrmv_NUN, dtrmv_NLU, dtrmv_NLN,
  dtrmv_TUU, dtrmv_TUN, dtrmv_TLU, dtrmv_TLN,
#නැතහොත්
  strmv_NUU, strmv_NUN, strmv_NLU, strmv_NLN,
  strmv_TUU, strmv_TUN, strmv_TLU, strmv_TLN,
#endif
};

#ifdef SMP
static පූර්ණ_සංඛ්‍යා (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, පූර්ණ_සංඛ්‍යා) = {
#ifdef XDOUBLE
  qtrmv_thread_NUU, qtrmv_thread_NUN, qtrmv_thread_NLU, qtrmv_thread_NLN,
  qtrmv_thread_TUU, qtrmv_thread_TUN, qtrmv_thread_TLU, qtrmv_thread_TLN,
#elif defined(DOUBLE)
  dtrmv_thread_NUU, dtrmv_thread_NUN, dtrmv_thread_NLU, dtrmv_thread_NLN,
  dtrmv_thread_TUU, dtrmv_thread_TUN, dtrmv_thread_TLU, dtrmv_thread_TLN,
#නැතහොත්
  strmv_thread_NUU, strmv_thread_NUN, strmv_thread_NLU, strmv_thread_NLN,
  strmv_thread_TUU, strmv_thread_TUN, strmv_thread_TLU, strmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANS, var_na *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;
  var_na diag_arg  = *DIAG;

  blasint n    = *N;
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

  නම් (incx == 0)          info =  8;
  නම් (lda  < MAX(1, n))   info =  6;
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
	   blasint n, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

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

    නම් (incx == 0)          info =  8;
    නම් (lda  < MAX(1, n))   info =  6;
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

    නම් (incx == 0)          info =  8;
    නම් (lda  < MAX(1, n))   info =  6;
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
  
  නම් (n < 50 ) nthreads = 1;
  නම් (nthreads > 2 && n < 500) nthreads = 2;
	
  නම් (nthreads == 1) {
#endif

    (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

#ifdef SMP
  } නැතහොත් {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n,  n * n);

  IDEBUG_END;

  ආපසු;
}
