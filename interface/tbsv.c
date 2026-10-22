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
#ꦭꦶꦪ
#define ERROR_NAME "STBSV "
#endif

static puur_naa_nka (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtbsv_NUU, qtbsv_NUN, qtbsv_NLU, qtbsv_NLN,
  qtbsv_TUU, qtbsv_TUN, qtbsv_TLU, qtbsv_TLN,
#elif defined(DOUBLE)
  dtbsv_NUU, dtbsv_NUN, dtbsv_NLU, dtbsv_NLN,
  dtbsv_TUU, dtbsv_TUN, dtbsv_TLU, dtbsv_TLN,
#ꦭꦶꦪ
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
  puur_naa_nka uplo;
  puur_naa_nka unit;
  puur_naa_nka trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ꦪꦼꦤ꧀ (trans_arg == 'N') trans = 0;
  ꦪꦼꦤ꧀ (trans_arg == 'T') trans = 1;
  ꦪꦼꦤ꧀ (trans_arg == 'R') trans = 0;
  ꦪꦼꦤ꧀ (trans_arg == 'C') trans = 1;

  ꦪꦼꦤ꧀ (diag_arg  == 'U') unit  = 0;
  ꦪꦼꦤ꧀ (diag_arg  == 'N') unit  = 1;

  ꦪꦼꦤ꧀ (uplo_arg  == 'U') uplo  = 0;
  ꦪꦼꦤ꧀ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ꦪꦼꦤ꧀ (incx == 0)          info =  9;
  ꦪꦼꦤ꧀ (lda < k + 1)        info =  7;
  ꦪꦼꦤ꧀ (k < 0)              info =  5;
  ꦪꦼꦤ꧀ (n < 0)              info =  4;
  ꦪꦼꦤ꧀ (unit  < 0)          info =  3;
  ꦪꦼꦤ꧀ (trans < 0)          info =  2;
  ꦪꦼꦤ꧀ (uplo  < 0)          info =  1;

  ꦪꦼꦤ꧀ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ꦮꦭꦶ;
  }

#ꦭꦶꦪ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  puur_naa_nka trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ꦪꦼꦤ꧀ (order == CblasColMajor) {
    ꦪꦼꦤ꧀ (Uplo == CblasUpper)         uplo  = 0;
    ꦪꦼꦤ꧀ (Uplo == CblasLower)         uplo  = 1;

    ꦪꦼꦤ꧀ (TransA == CblasNoTrans)     trans = 0;
    ꦪꦼꦤ꧀ (TransA == CblasTrans)       trans = 1;
    ꦪꦼꦤ꧀ (TransA == CblasConjNoTrans) trans = 0;
    ꦪꦼꦤ꧀ (TransA == CblasConjTrans)   trans = 1;

    ꦪꦼꦤ꧀ (Diag == CblasUnit)          unit  = 0;
    ꦪꦼꦤ꧀ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ꦪꦼꦤ꧀ (incx == 0)          info =  9;
    ꦪꦼꦤ꧀ (lda < k + 1)        info =  7;
    ꦪꦼꦤ꧀ (k < 0)              info =  5;
    ꦪꦼꦤ꧀ (n < 0)              info =  4;
    ꦪꦼꦤ꧀ (unit  < 0)          info =  3;
    ꦪꦼꦤ꧀ (trans < 0)          info =  2;
    ꦪꦼꦤ꧀ (uplo  < 0)          info =  1;
  }

  ꦪꦼꦤ꧀ (order == CblasRowMajor) {
    ꦪꦼꦤ꧀ (Uplo == CblasUpper)         uplo  = 1;
    ꦪꦼꦤ꧀ (Uplo == CblasLower)         uplo  = 0;

    ꦪꦼꦤ꧀ (TransA == CblasNoTrans)     trans = 1;
    ꦪꦼꦤ꧀ (TransA == CblasTrans)       trans = 0;
    ꦪꦼꦤ꧀ (TransA == CblasConjNoTrans) trans = 1;
    ꦪꦼꦤ꧀ (TransA == CblasConjTrans)   trans = 0;

    ꦪꦼꦤ꧀ (Diag == CblasUnit)          unit  = 0;
    ꦪꦼꦤ꧀ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ꦪꦼꦤ꧀ (incx == 0)          info =  9;
    ꦪꦼꦤ꧀ (lda < k + 1)        info =  7;
    ꦪꦼꦤ꧀ (k < 0)              info =  5;
    ꦪꦼꦤ꧀ (n < 0)              info =  4;
    ꦪꦼꦤ꧀ (unit  < 0)          info =  3;
    ꦪꦼꦤ꧀ (trans < 0)          info =  2;
    ꦪꦼꦤ꧀ (uplo  < 0)          info =  1;
  }

  ꦪꦼꦤ꧀ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ꦮꦭꦶ;
  }

#endif

  ꦪꦼꦤ꧀ (n == 0) ꦮꦭꦶ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ꦪꦼꦤ꧀ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ꦪꦼꦤ꧀ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ꦮꦭꦶ;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  ꦮꦭꦶ;
}
