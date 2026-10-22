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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <ctype.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XTPSV "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZTPSV "
#ଅନ୍ୟଥା
#define ERROR_NAME "CTPSV "
#endif

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, ଶୂନ୍ୟ *) = {
#ifdef XDOUBLE
  xtpsv_NUU, xtpsv_NUN, xtpsv_NLU, xtpsv_NLN,
  xtpsv_TUU, xtpsv_TUN, xtpsv_TLU, xtpsv_TLN,
  xtpsv_RUU, xtpsv_RUN, xtpsv_RLU, xtpsv_RLN,
  xtpsv_CUU, xtpsv_CUN, xtpsv_CLU, xtpsv_CLN,
#ନଚେତଯଦି defined(DOUBLE)
  ztpsv_NUU, ztpsv_NUN, ztpsv_NLU, ztpsv_NLN,
  ztpsv_TUU, ztpsv_TUN, ztpsv_TLU, ztpsv_TLN,
  ztpsv_RUU, ztpsv_RUN, ztpsv_RLU, ztpsv_RLN,
  ztpsv_CUU, ztpsv_CUN, ztpsv_CLU, ztpsv_CLN,
#ଅନ୍ୟଥା
  ctpsv_NUU, ctpsv_NUN, ctpsv_NLU, ctpsv_NLN,
  ctpsv_TUU, ctpsv_TUN, ctpsv_TLU, ctpsv_TLN,
  ctpsv_RUU, ctpsv_RUN, ctpsv_RLU, ctpsv_RLN,
  ctpsv_CUU, ctpsv_CUN, ctpsv_CLU, ctpsv_CLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, ଅକ୍ଷର *TRANS, ଅକ୍ଷର *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  ଅକ୍ଷର uplo_arg  = *UPLO;
  ଅକ୍ଷର trans_arg = *TRANS;
  ଅକ୍ଷର diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା unit;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'T') trans = 1;
  ଯଦି (trans_arg == 'R') trans = 2;
  ଯଦି (trans_arg == 'C') trans = 3;

  ଯଦି (diag_arg  == 'U') unit  = 0;
  ଯଦି (diag_arg  == 'N') unit  = 1;

  ଯଦି (uplo_arg  == 'U') uplo  = 0;
  ଯଦି (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ଯଦି (incx == 0)          info =  7;
  ଯଦି (n < 0)              info =  4;
  ଯଦି (unit  < 0)          info =  3;
  ଯଦି (trans < 0)          info =  2;
  ଯଦି (uplo  < 0)          info =  1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_UPLO Uplo,
	   ଗଣନା CBLAS_TRANSPOSE TransA, ଗଣନା CBLAS_DIAG Diag,
	   blasint n, ଶୂନ୍ୟ  *va, ଶୂନ୍ୟ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (Uplo == CblasUpper)         uplo  = 0;
    ଯଦି (Uplo == CblasLower)         uplo  = 1;

    ଯଦି (TransA == CblasNoTrans)     trans = 0;
    ଯଦି (TransA == CblasTrans)       trans = 1;
    ଯଦି (TransA == CblasConjNoTrans) trans = 2;
    ଯଦି (TransA == CblasConjTrans)   trans = 3;

    ଯଦି (Diag == CblasUnit)          unit  = 0;
    ଯଦି (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ଯଦି (incx == 0)          info =  7;
    ଯଦି (n < 0)              info =  4;
    ଯଦି (unit  < 0)          info =  3;
    ଯଦି (trans < 0)          info =  2;
    ଯଦି (uplo  < 0)          info =  1;
  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (Uplo == CblasUpper)         uplo  = 1;
    ଯଦି (Uplo == CblasLower)         uplo  = 0;

    ଯଦି (TransA == CblasNoTrans)     trans = 1;
    ଯଦି (TransA == CblasTrans)       trans = 0;
    ଯଦି (TransA == CblasConjNoTrans) trans = 3;
    ଯଦି (TransA == CblasConjTrans)   trans = 2;

    ଯଦି (Diag == CblasUnit)          unit  = 0;
    ଯଦି (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ଯଦି (incx == 0)          info =  7;
    ଯଦି (n < 0)              info =  4;
    ଯଦି (unit  < 0)          info =  3;
    ଯଦି (trans < 0)          info =  2;
    ଯଦି (uplo  < 0)          info =  1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#endif

  ଯଦି (n == 0) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ଯଦି (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
