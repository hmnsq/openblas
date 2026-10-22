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
#define ERROR_NAME "XTRMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTRMV "
#ᤀᤡ᤺ᤆᤡ
#define ERROR_NAME "CTRMV "
#endif

static puur_naa_nka (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xtrmv_NUU, xtrmv_NUN, xtrmv_NLU, xtrmv_NLN,
  xtrmv_TUU, xtrmv_TUN, xtrmv_TLU, xtrmv_TLN,
  xtrmv_RUU, xtrmv_RUN, xtrmv_RLU, xtrmv_RLN,
  xtrmv_CUU, xtrmv_CUN, xtrmv_CLU, xtrmv_CLN,
#elif defined(DOUBLE)
  ztrmv_NUU, ztrmv_NUN, ztrmv_NLU, ztrmv_NLN,
  ztrmv_TUU, ztrmv_TUN, ztrmv_TLU, ztrmv_TLN,
  ztrmv_RUU, ztrmv_RUN, ztrmv_RLU, ztrmv_RLN,
  ztrmv_CUU, ztrmv_CUN, ztrmv_CLU, ztrmv_CLN,
#ᤀᤡ᤺ᤆᤡ
  ctrmv_NUU, ctrmv_NUN, ctrmv_NLU, ctrmv_NLN,
  ctrmv_TUU, ctrmv_TUN, ctrmv_TLU, ctrmv_TLN,
  ctrmv_RUU, ctrmv_RUN, ctrmv_RLU, ctrmv_RLN,
  ctrmv_CUU, ctrmv_CUN, ctrmv_CLU, ctrmv_CLN,
#endif
};

#ifdef SMP
static puur_naa_nka (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  xtrmv_thread_NUU, xtrmv_thread_NUN, xtrmv_thread_NLU, xtrmv_thread_NLN,
  xtrmv_thread_TUU, xtrmv_thread_TUN, xtrmv_thread_TLU, xtrmv_thread_TLN,
  xtrmv_thread_RUU, xtrmv_thread_RUN, xtrmv_thread_RLU, xtrmv_thread_RLN,
  xtrmv_thread_CUU, xtrmv_thread_CUN, xtrmv_thread_CLU, xtrmv_thread_CLN,
#elif defined(DOUBLE)
  ztrmv_thread_NUU, ztrmv_thread_NUN, ztrmv_thread_NLU, ztrmv_thread_NLN,
  ztrmv_thread_TUU, ztrmv_thread_TUN, ztrmv_thread_TLU, ztrmv_thread_TLN,
  ztrmv_thread_RUU, ztrmv_thread_RUN, ztrmv_thread_RLU, ztrmv_thread_RLN,
  ztrmv_thread_CUU, ztrmv_thread_CUN, ztrmv_thread_CLU, ztrmv_thread_CLN,
#ᤀᤡ᤺ᤆᤡ
  ctrmv_thread_NUU, ctrmv_thread_NUN, ctrmv_thread_NLU, ctrmv_thread_NLN,
  ctrmv_thread_TUU, ctrmv_thread_TUN, ctrmv_thread_TLU, ctrmv_thread_TLN,
  ctrmv_thread_RUU, ctrmv_thread_RUN, ctrmv_thread_RLU, ctrmv_thread_RLN,
  ctrmv_thread_CUU, ctrmv_thread_CUN, ctrmv_thread_CLU, ctrmv_thread_CLN,
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
  puur_naa_nka uplo;
  puur_naa_nka unit;
  puur_naa_nka trans, buffer_size;
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

  ᤀᤡ᤺ (trans_arg == 'N') trans = 0;
  ᤀᤡ᤺ (trans_arg == 'T') trans = 1;
  ᤀᤡ᤺ (trans_arg == 'R') trans = 2;
  ᤀᤡ᤺ (trans_arg == 'C') trans = 3;

  ᤀᤡ᤺ (diag_arg  == 'U') unit  = 0;
  ᤀᤡ᤺ (diag_arg  == 'N') unit  = 1;

  ᤀᤡ᤺ (uplo_arg  == 'U') uplo  = 0;
  ᤀᤡ᤺ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᤀᤡ᤺ (incx == 0)          info =  8;
  ᤀᤡ᤺ (lda  < MAX(1, n))   info =  6;
  ᤀᤡ᤺ (n < 0)              info =  4;
  ᤀᤡ᤺ (unit  < 0)          info =  3;
  ᤀᤡ᤺ (trans < 0)          info =  2;
  ᤀᤡ᤺ (uplo  < 0)          info =  1;

  ᤀᤡ᤺ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#ᤀᤡ᤺ᤆᤡ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, void  *va, blasint lda, void  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  puur_naa_nka trans, uplo, unit, buffer_size;
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

  ᤀᤡ᤺ (order == CblasColMajor) {
    ᤀᤡ᤺ (Uplo == CblasUpper)         uplo  = 0;
    ᤀᤡ᤺ (Uplo == CblasLower)         uplo  = 1;

    ᤀᤡ᤺ (TransA == CblasNoTrans)     trans = 0;
    ᤀᤡ᤺ (TransA == CblasTrans)       trans = 1;
    ᤀᤡ᤺ (TransA == CblasConjNoTrans) trans = 2;
    ᤀᤡ᤺ (TransA == CblasConjTrans)   trans = 3;

    ᤀᤡ᤺ (Diag == CblasUnit)          unit  = 0;
    ᤀᤡ᤺ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ᤀᤡ᤺ (incx == 0)          info =  8;
    ᤀᤡ᤺ (lda  < MAX(1, n))   info =  6;
    ᤀᤡ᤺ (n < 0)              info =  4;
    ᤀᤡ᤺ (unit  < 0)          info =  3;
    ᤀᤡ᤺ (trans < 0)          info =  2;
    ᤀᤡ᤺ (uplo  < 0)          info =  1;
  }

  ᤀᤡ᤺ (order == CblasRowMajor) {
    ᤀᤡ᤺ (Uplo == CblasUpper)         uplo  = 1;
    ᤀᤡ᤺ (Uplo == CblasLower)         uplo  = 0;

    ᤀᤡ᤺ (TransA == CblasNoTrans)     trans = 1;
    ᤀᤡ᤺ (TransA == CblasTrans)       trans = 0;
    ᤀᤡ᤺ (TransA == CblasConjNoTrans) trans = 3;
    ᤀᤡ᤺ (TransA == CblasConjTrans)   trans = 2;

    ᤀᤡ᤺ (Diag == CblasUnit)          unit  = 0;
    ᤀᤡ᤺ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ᤀᤡ᤺ (incx == 0)          info =  8;
    ᤀᤡ᤺ (lda  < MAX(1, n))   info =  6;
    ᤀᤡ᤺ (n < 0)              info =  4;
    ᤀᤡ᤺ (unit  < 0)          info =  3;
    ᤀᤡ᤺ (trans < 0)          info =  2;
    ᤀᤡ᤺ (uplo  < 0)          info =  1;
  }

  ᤀᤡ᤺ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#endif

  ᤀᤡ᤺ (n == 0) ᤀᤡᤴ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᤀᤡ᤺ (incx < 0 ) x -= (n - 1) * incx * 2;

#ifdef SMP
  // Calibrated on a Xeon E5-2630
  ᤀᤡ᤺(1L * n * n > 36L * parimaa_na(FLOAT) * parimaa_na(FLOAT) * GEMM_MULTITHREAD_THRESHOLD) {
    nthreads = num_cpu_avail(2);
    ᤀᤡ᤺(nthreads > 2 && 1L * n * n < 64L * parimaa_na(FLOAT) * parimaa_na(FLOAT) * GEMM_MULTITHREAD_THRESHOLD)
      nthreads = 2;
  } ᤀᤡ᤺ᤆᤡ
      nthreads = 1;

  ᤀᤡ᤺(nthreads > 1) {
    buffer_size = n > 16 ? 0 : n * 4 + 40;
  }
  ᤀᤡ᤺ᤆᤡ
#endif
  {
    buffer_size = ((n - 1) / DTB_ENTRIES) * 2 * DTB_ENTRIES + 32 / parimaa_na(FLOAT);
    // It seems to be required for some K8 or Barcelona CPU
    buffer_size += 8;
    ᤀᤡ᤺(incx != 1)
      buffer_size += n * 2;
  }
  STACK_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  ᤀᤡ᤺ (nthreads == 1) {
#endif

  (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

#ifdef SMP
  } ᤀᤡ᤺ᤆᤡ {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer, nthreads);

  }
#endif

  STACK_FREE(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n,  n * n);

  IDEBUG_END;

  ᤀᤡᤴ;
}
