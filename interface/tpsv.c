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
#define ERROR_NAME "QTPSV "
#elif defined(DOUBLE)
#define ERROR_NAME "DTPSV "
#မဟုတ်ပါက
#define ERROR_NAME "STPSV "
#endif

static ကိန်းပြည့် (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtpsv_NUU, qtpsv_NUN, qtpsv_NLU, qtpsv_NLN,
  qtpsv_TUU, qtpsv_TUN, qtpsv_TLU, qtpsv_TLN,
#elif defined(DOUBLE)
  dtpsv_NUU, dtpsv_NUN, dtpsv_NLU, dtpsv_NLN,
  dtpsv_TUU, dtpsv_TUN, dtpsv_TLU, dtpsv_TLN,
#မဟုတ်ပါက
  stpsv_NUU, stpsv_NUN, stpsv_NLU, stpsv_NLN,
  stpsv_TUU, stpsv_TUN, stpsv_TLU, stpsv_TLN,
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
  ကိန်းပြည့် uplo;
  ကိန်းပြည့် unit;
  ကိန်းပြည့် trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  အကယ်၍ (trans_arg == 'N') trans = 0;
  အကယ်၍ (trans_arg == 'T') trans = 1;
  အကယ်၍ (trans_arg == 'R') trans = 0;
  အကယ်၍ (trans_arg == 'C') trans = 1;

  အကယ်၍ (diag_arg  == 'U') unit  = 0;
  အကယ်၍ (diag_arg  == 'N') unit  = 1;

  အကယ်၍ (uplo_arg  == 'U') uplo  = 0;
  အကယ်၍ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  အကယ်၍ (incx == 0)          info =  7;
  အကယ်၍ (n < 0)              info =  4;
  အကယ်၍ (unit  < 0)          info =  3;
  အကယ်၍ (trans < 0)          info =  2;
  အကယ်၍ (uplo  < 0)          info =  1;

  အကယ်၍ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#မဟုတ်ပါက

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {

  ကိန်းပြည့် trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  အကယ်၍ (order == CblasColMajor) {
    အကယ်၍ (Uplo == CblasUpper)         uplo  = 0;
    အကယ်၍ (Uplo == CblasLower)         uplo  = 1;

    အကယ်၍ (TransA == CblasNoTrans)     trans = 0;
    အကယ်၍ (TransA == CblasTrans)       trans = 1;
    အကယ်၍ (TransA == CblasConjNoTrans) trans = 0;
    အကယ်၍ (TransA == CblasConjTrans)   trans = 1;

    အကယ်၍ (Diag == CblasUnit)          unit  = 0;
    အကယ်၍ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    အကယ်၍ (incx == 0)          info =  7;
    အကယ်၍ (n < 0)              info =  4;
    အကယ်၍ (unit  < 0)          info =  3;
    အကယ်၍ (trans < 0)          info =  2;
    အကယ်၍ (uplo  < 0)          info =  1;
  }

  အကယ်၍ (order == CblasRowMajor) {
    အကယ်၍ (Uplo == CblasUpper)         uplo  = 1;
    အကယ်၍ (Uplo == CblasLower)         uplo  = 0;

    အကယ်၍ (TransA == CblasNoTrans)     trans = 1;
    အကယ်၍ (TransA == CblasTrans)       trans = 0;
    အကယ်၍ (TransA == CblasConjNoTrans) trans = 1;
    အကယ်၍ (TransA == CblasConjTrans)   trans = 0;

    အကယ်၍ (Diag == CblasUnit)          unit  = 0;
    အကယ်၍ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    အကယ်၍ (incx == 0)          info =  7;
    အကယ်၍ (n < 0)              info =  4;
    အကယ်၍ (unit  < 0)          info =  3;
    အကယ်၍ (trans < 0)          info =  2;
    အကယ်၍ (uplo  < 0)          info =  1;
  }

  အကယ်၍ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#endif

  အကယ်၍ (n == 0) ပြန်ပေး;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  အကယ်၍ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  အကယ်၍ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  ပြန်ပေး;
}
