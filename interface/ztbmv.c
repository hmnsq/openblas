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
#define ERROR_NAME "XTBMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTBMV "
#มิฉะนั้น
#define ERROR_NAME "CTBMV "
#endif

static จำนวนเต็ม (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  xtbmv_NUU, xtbmv_NUN, xtbmv_NLU, xtbmv_NLN,
  xtbmv_TUU, xtbmv_TUN, xtbmv_TLU, xtbmv_TLN,
  xtbmv_RUU, xtbmv_RUN, xtbmv_RLU, xtbmv_RLN,
  xtbmv_CUU, xtbmv_CUN, xtbmv_CLU, xtbmv_CLN,
#elif defined(DOUBLE)
  ztbmv_NUU, ztbmv_NUN, ztbmv_NLU, ztbmv_NLN,
  ztbmv_TUU, ztbmv_TUN, ztbmv_TLU, ztbmv_TLN,
  ztbmv_RUU, ztbmv_RUN, ztbmv_RLU, ztbmv_RLN,
  ztbmv_CUU, ztbmv_CUN, ztbmv_CLU, ztbmv_CLN,
#มิฉะนั้น
  ctbmv_NUU, ctbmv_NUN, ctbmv_NLU, ctbmv_NLN,
  ctbmv_TUU, ctbmv_TUN, ctbmv_TLU, ctbmv_TLN,
  ctbmv_RUU, ctbmv_RUN, ctbmv_RLU, ctbmv_RLN,
  ctbmv_CUU, ctbmv_CUN, ctbmv_CLU, ctbmv_CLN,
#endif
};

#ifdef SMP
static จำนวนเต็ม (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, จำนวนเต็ม) = {
#ifdef XDOUBLE
  xtbmv_thread_NUU, xtbmv_thread_NUN, xtbmv_thread_NLU, xtbmv_thread_NLN,
  xtbmv_thread_TUU, xtbmv_thread_TUN, xtbmv_thread_TLU, xtbmv_thread_TLN,
  xtbmv_thread_RUU, xtbmv_thread_RUN, xtbmv_thread_RLU, xtbmv_thread_RLN,
  xtbmv_thread_CUU, xtbmv_thread_CUN, xtbmv_thread_CLU, xtbmv_thread_CLN,
#elif defined(DOUBLE)
  ztbmv_thread_NUU, ztbmv_thread_NUN, ztbmv_thread_NLU, ztbmv_thread_NLN,
  ztbmv_thread_TUU, ztbmv_thread_TUN, ztbmv_thread_TLU, ztbmv_thread_TLN,
  ztbmv_thread_RUU, ztbmv_thread_RUN, ztbmv_thread_RLU, ztbmv_thread_RLN,
  ztbmv_thread_CUU, ztbmv_thread_CUN, ztbmv_thread_CLU, ztbmv_thread_CLN,
#มิฉะนั้น
  ctbmv_thread_NUU, ctbmv_thread_NUN, ctbmv_thread_NLU, ctbmv_thread_NLN,
  ctbmv_thread_TUU, ctbmv_thread_TUN, ctbmv_thread_TLU, ctbmv_thread_TLN,
  ctbmv_thread_RUU, ctbmv_thread_RUN, ctbmv_thread_RLU, ctbmv_thread_RLN,
  ctbmv_thread_CUU, ctbmv_thread_CUN, ctbmv_thread_CLU, ctbmv_thread_CLN,
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
  จำนวนเต็ม uplo;
  จำนวนเต็ม unit;
  จำนวนเต็ม trans;
  FLOAT *buffer;
#ifdef SMP
  จำนวนเต็ม nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ถ้า (trans_arg == 'N') trans = 0;
  ถ้า (trans_arg == 'T') trans = 1;
  ถ้า (trans_arg == 'R') trans = 2;
  ถ้า (trans_arg == 'C') trans = 3;

  ถ้า (diag_arg  == 'U') unit  = 0;
  ถ้า (diag_arg  == 'N') unit  = 1;

  ถ้า (uplo_arg  == 'U') uplo  = 0;
  ถ้า (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ถ้า (incx == 0)          info =  9;
  ถ้า (lda < k + 1)        info =  7;
  ถ้า (k < 0)              info =  5;
  ถ้า (n < 0)              info =  4;
  ถ้า (unit  < 0)          info =  3;
  ถ้า (trans < 0)          info =  2;
  ถ้า (uplo  < 0)          info =  1;

  ถ้า (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#มิฉะนั้น

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, blasint k, void  *va, blasint lda, void  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;

  จำนวนเต็ม trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  จำนวนเต็ม nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ถ้า (order == CblasColMajor) {
    ถ้า (Uplo == CblasUpper)         uplo  = 0;
    ถ้า (Uplo == CblasLower)         uplo  = 1;

    ถ้า (TransA == CblasNoTrans)     trans = 0;
    ถ้า (TransA == CblasTrans)       trans = 1;
    ถ้า (TransA == CblasConjNoTrans) trans = 2;
    ถ้า (TransA == CblasConjTrans)   trans = 3;

    ถ้า (Diag == CblasUnit)          unit  = 0;
    ถ้า (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ถ้า (incx == 0)          info =  9;
    ถ้า (lda < k + 1)        info =  7;
    ถ้า (k < 0)              info =  5;
    ถ้า (n < 0)              info =  4;
    ถ้า (unit  < 0)          info =  3;
    ถ้า (trans < 0)          info =  2;
    ถ้า (uplo  < 0)          info =  1;
  }

  ถ้า (order == CblasRowMajor) {
    ถ้า (Uplo == CblasUpper)         uplo  = 1;
    ถ้า (Uplo == CblasLower)         uplo  = 0;

    ถ้า (TransA == CblasNoTrans)     trans = 1;
    ถ้า (TransA == CblasTrans)       trans = 0;
    ถ้า (TransA == CblasConjNoTrans) trans = 3;
    ถ้า (TransA == CblasConjTrans)   trans = 2;

    ถ้า (Diag == CblasUnit)          unit  = 0;
    ถ้า (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ถ้า (incx == 0)          info =  9;
    ถ้า (lda < k + 1)        info =  7;
    ถ้า (k < 0)              info =  5;
    ถ้า (n < 0)              info =  4;
    ถ้า (unit  < 0)          info =  3;
    ถ้า (trans < 0)          info =  2;
    ถ้า (uplo  < 0)          info =  1;
  }

  ถ้า (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#endif

  ถ้า (n == 0) คืน;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ถ้า (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ถ้า (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ถ้า (n < 50) nthreads = 1;
  ถ้า (nthreads > 2 && n < 500) nthreads = 2;

  ถ้า (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } มิฉะนั้น {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  คืน;
}
