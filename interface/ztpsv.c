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
#define ERROR_NAME "XTPSV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTPSV "
#បើមិនដូច្នេះ
#define ERROR_NAME "CTPSV "
#endif

static ចំនួនគត់ (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  xtpsv_NUU, xtpsv_NUN, xtpsv_NLU, xtpsv_NLN,
  xtpsv_TUU, xtpsv_TUN, xtpsv_TLU, xtpsv_TLN,
  xtpsv_RUU, xtpsv_RUN, xtpsv_RLU, xtpsv_RLN,
  xtpsv_CUU, xtpsv_CUN, xtpsv_CLU, xtpsv_CLN,
#elif defined(DOUBLE)
  ztpsv_NUU, ztpsv_NUN, ztpsv_NLU, ztpsv_NLN,
  ztpsv_TUU, ztpsv_TUN, ztpsv_TLU, ztpsv_TLN,
  ztpsv_RUU, ztpsv_RUN, ztpsv_RLU, ztpsv_RLN,
  ztpsv_CUU, ztpsv_CUN, ztpsv_CLU, ztpsv_CLN,
#បើមិនដូច្នេះ
  ctpsv_NUU, ctpsv_NUN, ctpsv_NLU, ctpsv_NLN,
  ctpsv_TUU, ctpsv_TUN, ctpsv_TLU, ctpsv_TLN,
  ctpsv_RUU, ctpsv_RUN, ctpsv_RLU, ctpsv_RLN,
  ctpsv_CUU, ctpsv_CUN, ctpsv_CLU, ctpsv_CLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANS, var_na *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;
  var_na diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  ចំនួនគត់ uplo;
  ចំនួនគត់ unit;
  ចំនួនគត់ trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ប្រសិនបើ (trans_arg == 'N') trans = 0;
  ប្រសិនបើ (trans_arg == 'T') trans = 1;
  ប្រសិនបើ (trans_arg == 'R') trans = 2;
  ប្រសិនបើ (trans_arg == 'C') trans = 3;

  ប្រសិនបើ (diag_arg  == 'U') unit  = 0;
  ប្រសិនបើ (diag_arg  == 'N') unit  = 1;

  ប្រសិនបើ (uplo_arg  == 'U') uplo  = 0;
  ប្រសិនបើ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ប្រសិនបើ (incx == 0)          info =  7;
  ប្រសិនបើ (n < 0)              info =  4;
  ប្រសិនបើ (unit  < 0)          info =  3;
  ប្រសិនបើ (trans < 0)          info =  2;
  ប្រសិនបើ (uplo  < 0)          info =  1;

  ប្រសិនបើ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

#បើមិនដូច្នេះ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, void  *va, void  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  ចំនួនគត់ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ប្រសិនបើ (order == CblasColMajor) {
    ប្រសិនបើ (Uplo == CblasUpper)         uplo  = 0;
    ប្រសិនបើ (Uplo == CblasLower)         uplo  = 1;

    ប្រសិនបើ (TransA == CblasNoTrans)     trans = 0;
    ប្រសិនបើ (TransA == CblasTrans)       trans = 1;
    ប្រសិនបើ (TransA == CblasConjNoTrans) trans = 2;
    ប្រសិនបើ (TransA == CblasConjTrans)   trans = 3;

    ប្រសិនបើ (Diag == CblasUnit)          unit  = 0;
    ប្រសិនបើ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ប្រសិនបើ (incx == 0)          info =  7;
    ប្រសិនបើ (n < 0)              info =  4;
    ប្រសិនបើ (unit  < 0)          info =  3;
    ប្រសិនបើ (trans < 0)          info =  2;
    ប្រសិនបើ (uplo  < 0)          info =  1;
  }

  ប្រសិនបើ (order == CblasRowMajor) {
    ប្រសិនបើ (Uplo == CblasUpper)         uplo  = 1;
    ប្រសិនបើ (Uplo == CblasLower)         uplo  = 0;

    ប្រសិនបើ (TransA == CblasNoTrans)     trans = 1;
    ប្រសិនបើ (TransA == CblasTrans)       trans = 0;
    ប្រសិនបើ (TransA == CblasConjNoTrans) trans = 3;
    ប្រសិនបើ (TransA == CblasConjTrans)   trans = 2;

    ប្រសិនបើ (Diag == CblasUnit)          unit  = 0;
    ប្រសិនបើ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ប្រសិនបើ (incx == 0)          info =  7;
    ប្រសិនបើ (n < 0)              info =  4;
    ប្រសិនបើ (unit  < 0)          info =  3;
    ប្រសិនបើ (trans < 0)          info =  2;
    ប្រសិនបើ (uplo  < 0)          info =  1;
  }

  ប្រសិនបើ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

#endif

  ប្រសិនបើ (n == 0) ត្រឡប់;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ប្រសិនបើ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ប្រសិនបើ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ត្រឡប់;
}
