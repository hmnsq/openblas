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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTBSV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DTBSV "
#אחרת
#define ERROR_NAME "STBSV "
#endif

ساكن صحيح (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
  qtbsv_NUU, qtbsv_NUN, qtbsv_NLU, qtbsv_NLN,
  qtbsv_TUU, qtbsv_TUN, qtbsv_TLU, qtbsv_TLN,
#وإلاإذا defined(DOUBLE)
  dtbsv_NUU, dtbsv_NUN, dtbsv_NLU, dtbsv_NLN,
  dtbsv_TUU, dtbsv_TUN, dtbsv_TLU, dtbsv_TLN,
#אחרת
  stbsv_NUU, stbsv_NUN, stbsv_NLU, stbsv_NLN,
  stbsv_TUU, stbsv_TUN, stbsv_TLU, stbsv_TLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  אם (trans_arg == 'N') trans = 0;
  אם (trans_arg == 'T') trans = 1;
  אם (trans_arg == 'R') trans = 0;
  אם (trans_arg == 'C') trans = 1;

  אם (diag_arg  == 'U') unit  = 0;
  אם (diag_arg  == 'N') unit  = 1;

  אם (uplo_arg  == 'U') uplo  = 0;
  אם (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  אם (incx == 0)          info =  9;
  אם (lda < k + 1)        info =  7;
  אם (k < 0)              info =  5;
  אם (n < 0)              info =  4;
  אם (unit  < 0)          info =  3;
  אם (trans < 0)          info =  2;
  אם (uplo  < 0)          info =  1;

  אם (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#אחרת

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  صحيح trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  אם (order == CblasColMajor) {
    אם (Uplo == CblasUpper)         uplo  = 0;
    אם (Uplo == CblasLower)         uplo  = 1;

    אם (TransA == CblasNoTrans)     trans = 0;
    אם (TransA == CblasTrans)       trans = 1;
    אם (TransA == CblasConjNoTrans) trans = 0;
    אם (TransA == CblasConjTrans)   trans = 1;

    אם (Diag == CblasUnit)          unit  = 0;
    אם (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    אם (incx == 0)          info =  9;
    אם (lda < k + 1)        info =  7;
    אם (k < 0)              info =  5;
    אם (n < 0)              info =  4;
    אם (unit  < 0)          info =  3;
    אם (trans < 0)          info =  2;
    אם (uplo  < 0)          info =  1;
  }

  אם (order == CblasRowMajor) {
    אם (Uplo == CblasUpper)         uplo  = 1;
    אם (Uplo == CblasLower)         uplo  = 0;

    אם (TransA == CblasNoTrans)     trans = 1;
    אם (TransA == CblasTrans)       trans = 0;
    אם (TransA == CblasConjNoTrans) trans = 1;
    אם (TransA == CblasConjTrans)   trans = 0;

    אם (Diag == CblasUnit)          unit  = 0;
    אם (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    אם (incx == 0)          info =  9;
    אם (lda < k + 1)        info =  7;
    אם (k < 0)              info =  5;
    אם (n < 0)              info =  4;
    אם (unit  < 0)          info =  3;
    אם (trans < 0)          info =  2;
    אם (uplo  < 0)          info =  1;
  }

  אם (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#endif

  אם (n == 0) החזר;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  אם (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  אם (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  החזר;
}
