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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTPSV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DTPSV "
#অন্যথা
#define ERROR_NAME "STPSV "
#endif

স্থির পূর্ণাঙ্ক (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
  qtpsv_NUU, qtpsv_NUN, qtpsv_NLU, qtpsv_NLN,
  qtpsv_TUU, qtpsv_TUN, qtpsv_TLU, qtpsv_TLN,
#নোচেদ্যদি defined(DOUBLE)
  dtpsv_NUU, dtpsv_NUN, dtpsv_NLU, dtpsv_NLN,
  dtpsv_TUU, dtpsv_TUN, dtpsv_TLU, dtpsv_TLN,
#অন্যথা
  stpsv_NUU, stpsv_NUN, stpsv_NLU, stpsv_NLN,
  stpsv_TUU, stpsv_TUN, stpsv_TLU, stpsv_TLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(বর্ণ *UPLO, বর্ণ *TRANS, বর্ণ *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  বর্ণ uplo_arg  = *UPLO;
  বর্ণ trans_arg = *TRANS;
  বর্ণ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  পূর্ণাঙ্ক uplo;
  পূর্ণাঙ্ক unit;
  পূর্ণাঙ্ক trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  যদি (trans_arg == 'N') trans = 0;
  যদি (trans_arg == 'T') trans = 1;
  যদি (trans_arg == 'R') trans = 0;
  যদি (trans_arg == 'C') trans = 1;

  যদি (diag_arg  == 'U') unit  = 0;
  যদি (diag_arg  == 'N') unit  = 1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  যদি (incx == 0)          info =  7;
  যদি (n < 0)              info =  4;
  যদি (unit  < 0)          info =  3;
  যদি (trans < 0)          info =  2;
  যদি (uplo  < 0)          info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo,
	   গণনা CBLAS_TRANSPOSE TransA, গণনা CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {

  পূর্ণাঙ্ক trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (Uplo == CblasUpper)         uplo  = 0;
    যদি (Uplo == CblasLower)         uplo  = 1;

    যদি (TransA == CblasNoTrans)     trans = 0;
    যদি (TransA == CblasTrans)       trans = 1;
    যদি (TransA == CblasConjNoTrans) trans = 0;
    যদি (TransA == CblasConjTrans)   trans = 1;

    যদি (Diag == CblasUnit)          unit  = 0;
    যদি (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    যদি (incx == 0)          info =  7;
    যদি (n < 0)              info =  4;
    যদি (unit  < 0)          info =  3;
    যদি (trans < 0)          info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (order == CblasRowMajor) {
    যদি (Uplo == CblasUpper)         uplo  = 1;
    যদি (Uplo == CblasLower)         uplo  = 0;

    যদি (TransA == CblasNoTrans)     trans = 1;
    যদি (TransA == CblasTrans)       trans = 0;
    যদি (TransA == CblasConjNoTrans) trans = 1;
    যদি (TransA == CblasConjTrans)   trans = 0;

    যদি (Diag == CblasUnit)          unit  = 0;
    যদি (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    যদি (incx == 0)          info =  7;
    যদি (n < 0)              info =  4;
    যদি (unit  < 0)          info =  3;
    যদি (trans < 0)          info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#endif

  যদি (n == 0) প্রত্যাগম;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  প্রত্যাগম;
}
