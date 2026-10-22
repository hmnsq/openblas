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
#define ERROR_NAME "XTPSV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZTPSV "
#وإلا
#define ERROR_NAME "CTPSV "
#endif

ساكن صحيح (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
  xtpsv_NUU, xtpsv_NUN, xtpsv_NLU, xtpsv_NLN,
  xtpsv_TUU, xtpsv_TUN, xtpsv_TLU, xtpsv_TLN,
  xtpsv_RUU, xtpsv_RUN, xtpsv_RLU, xtpsv_RLN,
  xtpsv_CUU, xtpsv_CUN, xtpsv_CLU, xtpsv_CLN,
#وإلاإذا defined(DOUBLE)
  ztpsv_NUU, ztpsv_NUN, ztpsv_NLU, ztpsv_NLN,
  ztpsv_TUU, ztpsv_TUN, ztpsv_TLU, ztpsv_TLN,
  ztpsv_RUU, ztpsv_RUN, ztpsv_RLU, ztpsv_RLN,
  ztpsv_CUU, ztpsv_CUN, ztpsv_CLU, ztpsv_CLN,
#وإلا
  ctpsv_NUU, ctpsv_NUN, ctpsv_NLU, ctpsv_NLN,
  ctpsv_TUU, ctpsv_TUN, ctpsv_TLU, ctpsv_TLN,
  ctpsv_RUU, ctpsv_RUN, ctpsv_RLU, ctpsv_RLN,
  ctpsv_CUU, ctpsv_CUN, ctpsv_CLU, ctpsv_CLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  FLOAT *buffer;

  اطبع_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (diag_arg  == 'U') unit  = 0;
  إذا (diag_arg  == 'N') unit  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (incx == 0)          info =  7;
  إذا (n < 0)              info =  4;
  إذا (unit  < 0)          info =  3;
  إذا (trans < 0)          info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_DIAG Diag,
	   blasint n, فراغ  *va, فراغ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  صحيح trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  اطبع_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  إذا (order == CblasColMajor) {
    إذا (Uplo == CblasUpper)         uplo  = 0;
    إذا (Uplo == CblasLower)         uplo  = 1;

    إذا (TransA == CblasNoTrans)     trans = 0;
    إذا (TransA == CblasTrans)       trans = 1;
    إذا (TransA == CblasConjNoTrans) trans = 2;
    إذا (TransA == CblasConjTrans)   trans = 3;

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    إذا (incx == 0)          info =  7;
    إذا (n < 0)              info =  4;
    إذا (unit  < 0)          info =  3;
    إذا (trans < 0)          info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (order == CblasRowMajor) {
    إذا (Uplo == CblasUpper)         uplo  = 1;
    إذا (Uplo == CblasLower)         uplo  = 0;

    إذا (TransA == CblasNoTrans)     trans = 1;
    إذا (TransA == CblasTrans)       trans = 0;
    إذا (TransA == CblasConjNoTrans) trans = 3;
    إذا (TransA == CblasConjTrans)   trans = 2;

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    إذا (incx == 0)          info =  7;
    إذا (n < 0)              info =  4;
    إذا (unit  < 0)          info =  3;
    إذا (trans < 0)          info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#endif

  إذا (n == 0) أرجع;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  إذا (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  أرجع;
}
