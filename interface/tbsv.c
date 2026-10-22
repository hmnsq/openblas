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
#define ERROR_NAME "QTBSV "
#elif defined(DOUBLE)
#define ERROR_NAME "DTBSV "
#มิฉะนั้น
#define ERROR_NAME "STBSV "
#endif

static จำนวนเต็ม (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtbsv_NUU, qtbsv_NUN, qtbsv_NLU, qtbsv_NLN,
  qtbsv_TUU, qtbsv_TUN, qtbsv_TLU, qtbsv_TLN,
#elif defined(DOUBLE)
  dtbsv_NUU, dtbsv_NUN, dtbsv_NLU, dtbsv_NLN,
  dtbsv_TUU, dtbsv_TUN, dtbsv_TLU, dtbsv_TLN,
#มิฉะนั้น
  stbsv_NUU, stbsv_NUN, stbsv_NLU, stbsv_NLN,
  stbsv_TUU, stbsv_TUN, stbsv_TLU, stbsv_TLN,
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
  ถ้า (trans_arg == 'R') trans = 0;
  ถ้า (trans_arg == 'C') trans = 1;

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
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

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
    ถ้า (TransA == CblasConjNoTrans) trans = 0;
    ถ้า (TransA == CblasConjTrans)   trans = 1;

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
    ถ้า (TransA == CblasConjNoTrans) trans = 1;
    ถ้า (TransA == CblasConjTrans)   trans = 0;

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

  ถ้า (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ถ้า (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  คืน;
}
