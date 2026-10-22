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
#ပေႃးဢမ်ႇ
#define ERROR_NAME "CTPSV "
#endif

static puur_naa_nka (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, void *) = {
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
#ပေႃးဢမ်ႇ
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

  ပေႃး (trans_arg == 'N') trans = 0;
  ပေႃး (trans_arg == 'T') trans = 1;
  ပေႃး (trans_arg == 'R') trans = 2;
  ပေႃး (trans_arg == 'C') trans = 3;

  ပေႃး (diag_arg  == 'U') unit  = 0;
  ပေႃး (diag_arg  == 'N') unit  = 1;

  ပေႃး (uplo_arg  == 'U') uplo  = 0;
  ပေႃး (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ပေႃး (incx == 0)          info =  7;
  ပေႃး (n < 0)              info =  4;
  ပေႃး (unit  < 0)          info =  3;
  ပေႃး (trans < 0)          info =  2;
  ပေႃး (uplo  < 0)          info =  1;

  ပေႃး (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

#ပေႃးဢမ်ႇ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, void  *va, void  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  puur_naa_nka trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ပေႃး (order == CblasColMajor) {
    ပေႃး (Uplo == CblasUpper)         uplo  = 0;
    ပေႃး (Uplo == CblasLower)         uplo  = 1;

    ပေႃး (TransA == CblasNoTrans)     trans = 0;
    ပေႃး (TransA == CblasTrans)       trans = 1;
    ပေႃး (TransA == CblasConjNoTrans) trans = 2;
    ပေႃး (TransA == CblasConjTrans)   trans = 3;

    ပေႃး (Diag == CblasUnit)          unit  = 0;
    ပေႃး (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ပေႃး (incx == 0)          info =  7;
    ပေႃး (n < 0)              info =  4;
    ပေႃး (unit  < 0)          info =  3;
    ပေႃး (trans < 0)          info =  2;
    ပေႃး (uplo  < 0)          info =  1;
  }

  ပေႃး (order == CblasRowMajor) {
    ပေႃး (Uplo == CblasUpper)         uplo  = 1;
    ပေႃး (Uplo == CblasLower)         uplo  = 0;

    ပေႃး (TransA == CblasNoTrans)     trans = 1;
    ပေႃး (TransA == CblasTrans)       trans = 0;
    ပေႃး (TransA == CblasConjNoTrans) trans = 3;
    ပေႃး (TransA == CblasConjTrans)   trans = 2;

    ပေႃး (Diag == CblasUnit)          unit  = 0;
    ပေႃး (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ပေႃး (incx == 0)          info =  7;
    ပေႃး (n < 0)              info =  4;
    ပေႃး (unit  < 0)          info =  3;
    ပေႃး (trans < 0)          info =  2;
    ပေႃး (uplo  < 0)          info =  1;
  }

  ပေႃး (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

#endif

  ပေႃး (n == 0) ပွၵ်ႈ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ပေႃး (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ပေႃး (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ပွၵ်ႈ;
}
