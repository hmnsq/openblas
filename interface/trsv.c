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
#define ERROR_NAME "QTRSV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DTRSV "
#وإلا
#define ERROR_NAME "STRSV "
#endif

ساكن صحيح (*trsv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
  qtrsv_NUU, qtrsv_NUN, qtrsv_NLU, qtrsv_NLN,
  qtrsv_TUU, qtrsv_TUN, qtrsv_TLU, qtrsv_TLN,
#وإلاإذا defined(DOUBLE)
  dtrsv_NUU, dtrsv_NUN, dtrsv_NLU, dtrsv_NLN,
  dtrsv_TUU, dtrsv_TUN, dtrsv_TLU, dtrsv_TLN,
#وإلا
  strsv_NUU, strsv_NUN, strsv_NLU, strsv_NLN,
  strsv_TUU, strsv_TUN, strsv_TLU, strsv_TLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
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

  𐤀𐤌 (trans_arg == 'N') trans = 0;
  𐤀𐤌 (trans_arg == 'T') trans = 1;
  𐤀𐤌 (trans_arg == 'R') trans = 0;
  𐤀𐤌 (trans_arg == 'C') trans = 1;

  𐤀𐤌 (diag_arg  == 'U') unit  = 0;
  𐤀𐤌 (diag_arg  == 'N') unit  = 1;

  𐤀𐤌 (uplo_arg  == 'U') uplo  = 0;
  𐤀𐤌 (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  𐤀𐤌 (incx == 0)          info =  8;
  𐤀𐤌 (lda  < MAX(1, n))   info =  6;
  𐤀𐤌 (n < 0)              info =  4;
  𐤀𐤌 (unit  < 0)          info =  3;
  𐤀𐤌 (trans < 0)          info =  2;
  𐤀𐤌 (uplo  < 0)          info =  1;

  𐤀𐤌 (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  صحيح trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  𐤀𐤌 (order == CblasColMajor) {
    𐤀𐤌 (Uplo == CblasUpper)         uplo  = 0;
    𐤀𐤌 (Uplo == CblasLower)         uplo  = 1;

    𐤀𐤌 (TransA == CblasNoTrans)     trans = 0;
    𐤀𐤌 (TransA == CblasTrans)       trans = 1;
    𐤀𐤌 (TransA == CblasConjNoTrans) trans = 0;
    𐤀𐤌 (TransA == CblasConjTrans)   trans = 1;

    𐤀𐤌 (Diag == CblasUnit)          unit  = 0;
    𐤀𐤌 (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    𐤀𐤌 (incx == 0)          info =  8;
    𐤀𐤌 (lda  < MAX(1, n))   info =  6;
    𐤀𐤌 (n < 0)              info =  4;
    𐤀𐤌 (unit  < 0)          info =  3;
    𐤀𐤌 (trans < 0)          info =  2;
    𐤀𐤌 (uplo  < 0)          info =  1;
  }

  𐤀𐤌 (order == CblasRowMajor) {
    𐤀𐤌 (Uplo == CblasUpper)         uplo  = 1;
    𐤀𐤌 (Uplo == CblasLower)         uplo  = 0;

    𐤀𐤌 (TransA == CblasNoTrans)     trans = 1;
    𐤀𐤌 (TransA == CblasTrans)       trans = 0;
    𐤀𐤌 (TransA == CblasConjNoTrans) trans = 1;
    𐤀𐤌 (TransA == CblasConjTrans)   trans = 0;

    𐤀𐤌 (Diag == CblasUnit)          unit  = 0;
    𐤀𐤌 (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    𐤀𐤌 (incx == 0)          info =  8;
    𐤀𐤌 (lda  < MAX(1, n))   info =  6;
    𐤀𐤌 (n < 0)              info =  4;
    𐤀𐤌 (unit  < 0)          info =  3;
    𐤀𐤌 (trans < 0)          info =  2;
    𐤀𐤌 (uplo  < 0)          info =  1;
  }

  𐤀𐤌 (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#endif

  𐤀𐤌 (n == 0) أرجع;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  𐤀𐤌 (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  𐤀𐤌 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  (trsv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n,  n * n);

  IDEBUG_END;

  أرجع;
}
