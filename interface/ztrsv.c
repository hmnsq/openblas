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
#define ERROR_NAME "XTRSV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTRSV "
#ᮜᮤᮚ
#define ERROR_NAME "CTRSV "
#endif

static puur_naa_nka (*trsv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  xtrsv_NUU, xtrsv_NUN, xtrsv_NLU, xtrsv_NLN,
  xtrsv_TUU, xtrsv_TUN, xtrsv_TLU, xtrsv_TLN,
  xtrsv_RUU, xtrsv_RUN, xtrsv_RLU, xtrsv_RLN,
  xtrsv_CUU, xtrsv_CUN, xtrsv_CLU, xtrsv_CLN,
#elif defined(DOUBLE)
  ztrsv_NUU, ztrsv_NUN, ztrsv_NLU, ztrsv_NLN,
  ztrsv_TUU, ztrsv_TUN, ztrsv_TLU, ztrsv_TLN,
  ztrsv_RUU, ztrsv_RUN, ztrsv_RLU, ztrsv_RLN,
  ztrsv_CUU, ztrsv_CUN, ztrsv_CLU, ztrsv_CLN,
#ᮜᮤᮚ
  ctrsv_NUU, ctrsv_NUN, ctrsv_NLU, ctrsv_NLN,
  ctrsv_TUU, ctrsv_TUN, ctrsv_TLU, ctrsv_TLN,
  ctrsv_RUU, ctrsv_RUN, ctrsv_RLU, ctrsv_RLN,
  ctrsv_CUU, ctrsv_CUN, ctrsv_CLU, ctrsv_CLN,
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

  ᮚᮦᮔ᮪ (trans_arg == 'N') trans = 0;
  ᮚᮦᮔ᮪ (trans_arg == 'T') trans = 1;
  ᮚᮦᮔ᮪ (trans_arg == 'R') trans = 2;
  ᮚᮦᮔ᮪ (trans_arg == 'C') trans = 3;

  ᮚᮦᮔ᮪ (diag_arg  == 'U') unit  = 0;
  ᮚᮦᮔ᮪ (diag_arg  == 'N') unit  = 1;

  ᮚᮦᮔ᮪ (uplo_arg  == 'U') uplo  = 0;
  ᮚᮦᮔ᮪ (uplo_arg  == 'L') uplo  = 1;


  info = 0;

  ᮚᮦᮔ᮪ (incx == 0)          info =  8;
  ᮚᮦᮔ᮪ (lda  < MAX(1, n))   info =  6;
  ᮚᮦᮔ᮪ (n < 0)              info =  4;
  ᮚᮦᮔ᮪ (unit  < 0)          info =  3;
  ᮚᮦᮔ᮪ (trans < 0)          info =  2;
  ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;

  ᮚᮦᮔ᮪ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }


#ᮜᮤᮚ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, void  *va, blasint lda, void  *vx, blasint incx) {

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

  ᮚᮦᮔ᮪ (order == CblasColMajor) {
    ᮚᮦᮔ᮪ (Uplo == CblasUpper)         uplo  = 0;
    ᮚᮦᮔ᮪ (Uplo == CblasLower)         uplo  = 1;

    ᮚᮦᮔ᮪ (TransA == CblasNoTrans)     trans = 0;
    ᮚᮦᮔ᮪ (TransA == CblasTrans)       trans = 1;
    ᮚᮦᮔ᮪ (TransA == CblasConjNoTrans) trans = 2;
    ᮚᮦᮔ᮪ (TransA == CblasConjTrans)   trans = 3;

    ᮚᮦᮔ᮪ (Diag == CblasUnit)          unit  = 0;
    ᮚᮦᮔ᮪ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ᮚᮦᮔ᮪ (incx == 0)          info =  8;
    ᮚᮦᮔ᮪ (lda  < MAX(1, n))   info =  6;
    ᮚᮦᮔ᮪ (n < 0)              info =  4;
    ᮚᮦᮔ᮪ (unit  < 0)          info =  3;
    ᮚᮦᮔ᮪ (trans < 0)          info =  2;
    ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;
  }

  ᮚᮦᮔ᮪ (order == CblasRowMajor) {
    ᮚᮦᮔ᮪ (Uplo == CblasUpper)         uplo  = 1;
    ᮚᮦᮔ᮪ (Uplo == CblasLower)         uplo  = 0;

    ᮚᮦᮔ᮪ (TransA == CblasNoTrans)     trans = 1;
    ᮚᮦᮔ᮪ (TransA == CblasTrans)       trans = 0;
    ᮚᮦᮔ᮪ (TransA == CblasConjNoTrans) trans = 3;
    ᮚᮦᮔ᮪ (TransA == CblasConjTrans)   trans = 2;

    ᮚᮦᮔ᮪ (Diag == CblasUnit)          unit  = 0;
    ᮚᮦᮔ᮪ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ᮚᮦᮔ᮪ (incx == 0)          info =  8;
    ᮚᮦᮔ᮪ (lda  < MAX(1, n))   info =  6;
    ᮚᮦᮔ᮪ (n < 0)              info =  4;
    ᮚᮦᮔ᮪ (unit  < 0)          info =  3;
    ᮚᮦᮔ᮪ (trans < 0)          info =  2;
    ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;
  }

  ᮚᮦᮔ᮪ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#endif

  ᮚᮦᮔ᮪ (n == 0) ᮝᮜᮤ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᮚᮦᮔ᮪ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᮚᮦᮔ᮪ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

  (trsv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n,  n * n);

  IDEBUG_END;

  ᮝᮜᮤ;
}
