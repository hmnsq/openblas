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
#define ERROR_NAME "XTBSV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTBSV "
#มิฉะนั้น
#define ERROR_NAME "CTBSV "
#endif

static จำนวนเต็ม (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  xtbsv_NUU, xtbsv_NUN, xtbsv_NLU, xtbsv_NLN,
  xtbsv_TUU, xtbsv_TUN, xtbsv_TLU, xtbsv_TLN,
  xtbsv_RUU, xtbsv_RUN, xtbsv_RLU, xtbsv_RLN,
  xtbsv_CUU, xtbsv_CUN, xtbsv_CLU, xtbsv_CLN,
#elif defined(DOUBLE)
  ztbsv_NUU, ztbsv_NUN, ztbsv_NLU, ztbsv_NLN,
  ztbsv_TUU, ztbsv_TUN, ztbsv_TLU, ztbsv_TLN,
  ztbsv_RUU, ztbsv_RUN, ztbsv_RLU, ztbsv_RLN,
  ztbsv_CUU, ztbsv_CUN, ztbsv_CLU, ztbsv_CLN,
#มิฉะนั้น
  ctbsv_NUU, ctbsv_NUN, ctbsv_NLU, ctbsv_NLN,
  ctbsv_TUU, ctbsv_TUN, ctbsv_TLU, ctbsv_TLN,
  ctbsv_RUU, ctbsv_RUN, ctbsv_RLU, ctbsv_RLN,
  ctbsv_CUU, ctbsv_CUN, ctbsv_CLU, ctbsv_CLN,
#endif
};

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

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  คืน;
}
