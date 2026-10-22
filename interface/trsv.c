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
#define ERROR_NAME "QTRSV "
#elif defined(DOUBLE)
#define ERROR_NAME "DTRSV "
#དེ་མིན
#define ERROR_NAME "STRSV "
#endif

static ཧྲིལ (*trsv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtrsv_NUU, qtrsv_NUN, qtrsv_NLU, qtrsv_NLN,
  qtrsv_TUU, qtrsv_TUN, qtrsv_TLU, qtrsv_TLN,
#elif defined(DOUBLE)
  dtrsv_NUU, dtrsv_NUN, dtrsv_NLU, dtrsv_NLN,
  dtrsv_TUU, dtrsv_TUN, dtrsv_TLU, dtrsv_TLN,
#དེ་མིན
  strsv_NUU, strsv_NUN, strsv_NLU, strsv_NLN,
  strsv_TUU, strsv_TUN, strsv_TLU, strsv_TLN,
#endif
};

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
  ཧྲིལ uplo;
  ཧྲིལ unit;
  ཧྲིལ trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  གལ་སྲིད (trans_arg == 'N') trans = 0;
  གལ་སྲིད (trans_arg == 'T') trans = 1;
  གལ་སྲིད (trans_arg == 'R') trans = 0;
  གལ་སྲིད (trans_arg == 'C') trans = 1;

  གལ་སྲིད (diag_arg  == 'U') unit  = 0;
  གལ་སྲིད (diag_arg  == 'N') unit  = 1;

  གལ་སྲིད (uplo_arg  == 'U') uplo  = 0;
  གལ་སྲིད (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  གལ་སྲིད (incx == 0)          info =  8;
  གལ་སྲིད (lda  < MAX(1, n))   info =  6;
  གལ་སྲིད (n < 0)              info =  4;
  གལ་སྲིད (unit  < 0)          info =  3;
  གལ་སྲིད (trans < 0)          info =  2;
  གལ་སྲིད (uplo  < 0)          info =  1;

  གལ་སྲིད (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ལོག;
  }

#དེ་མིན

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  ཧྲིལ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  གལ་སྲིད (order == CblasColMajor) {
    གལ་སྲིད (Uplo == CblasUpper)         uplo  = 0;
    གལ་སྲིད (Uplo == CblasLower)         uplo  = 1;

    གལ་སྲིད (TransA == CblasNoTrans)     trans = 0;
    གལ་སྲིད (TransA == CblasTrans)       trans = 1;
    གལ་སྲིད (TransA == CblasConjNoTrans) trans = 0;
    གལ་སྲིད (TransA == CblasConjTrans)   trans = 1;

    གལ་སྲིད (Diag == CblasUnit)          unit  = 0;
    གལ་སྲིད (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    གལ་སྲིད (incx == 0)          info =  8;
    གལ་སྲིད (lda  < MAX(1, n))   info =  6;
    གལ་སྲིད (n < 0)              info =  4;
    གལ་སྲིད (unit  < 0)          info =  3;
    གལ་སྲིད (trans < 0)          info =  2;
    གལ་སྲིད (uplo  < 0)          info =  1;
  }

  གལ་སྲིད (order == CblasRowMajor) {
    གལ་སྲིད (Uplo == CblasUpper)         uplo  = 1;
    གལ་སྲིད (Uplo == CblasLower)         uplo  = 0;

    གལ་སྲིད (TransA == CblasNoTrans)     trans = 1;
    གལ་སྲིད (TransA == CblasTrans)       trans = 0;
    གལ་སྲིད (TransA == CblasConjNoTrans) trans = 1;
    གལ་སྲིད (TransA == CblasConjTrans)   trans = 0;

    གལ་སྲིད (Diag == CblasUnit)          unit  = 0;
    གལ་སྲིད (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    གལ་སྲིད (incx == 0)          info =  8;
    གལ་སྲིད (lda  < MAX(1, n))   info =  6;
    གལ་སྲིད (n < 0)              info =  4;
    གལ་སྲིད (unit  < 0)          info =  3;
    གལ་སྲིད (trans < 0)          info =  2;
    གལ་སྲིད (uplo  < 0)          info =  1;
  }

  གལ་སྲིད (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ལོག;
  }

#endif

  གལ་སྲིད (n == 0) ལོག;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  གལ་སྲིད (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  གལ་སྲིད (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ལོག;
  }

  (trsv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n,  n * n);

  IDEBUG_END;

  ལོག;
}
