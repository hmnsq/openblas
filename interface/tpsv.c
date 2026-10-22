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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <ctype.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTPSV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DTPSV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "STPSV "
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  qtpsv_NUU, qtpsv_NUN, qtpsv_NLU, qtpsv_NLN,
  qtpsv_TUU, qtpsv_TUN, qtpsv_TLU, qtpsv_TLN,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dtpsv_NUU, dtpsv_NUN, dtpsv_NLU, dtpsv_NLN,
  dtpsv_TUU, dtpsv_TUN, dtpsv_TLU, dtpsv_TLN,
#ಅನ್ಯಥಾ
  stpsv_NUU, stpsv_NUN, stpsv_NLU, stpsv_NLN,
  stpsv_TUU, stpsv_TUN, stpsv_TLU, stpsv_TLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, ವರ್ಣ *TRANS, ವರ್ಣ *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  ವರ್ಣ uplo_arg  = *UPLO;
  ವರ್ಣ trans_arg = *TRANS;
  ವರ್ಣ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  ಪೂರ್ಣಾಙ್ಕ unit;
  ಪೂರ್ಣಾಙ್ಕ trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
  ಯದಿ (trans_arg == 'R') trans = 0;
  ಯದಿ (trans_arg == 'C') trans = 1;

  ಯದಿ (diag_arg  == 'U') unit  = 0;
  ಯದಿ (diag_arg  == 'N') unit  = 1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ಯದಿ (incx == 0)          info =  7;
  ಯದಿ (n < 0)              info =  4;
  ಯದಿ (unit  < 0)          info =  3;
  ಯದಿ (trans < 0)          info =  2;
  ಯದಿ (uplo  < 0)          info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_UPLO Uplo,
	   ಗಣನಾ CBLAS_TRANSPOSE TransA, ಗಣನಾ CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {

  ಪೂರ್ಣಾಙ್ಕ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 0;
    ಯದಿ (Uplo == CblasLower)         uplo  = 1;

    ಯದಿ (TransA == CblasNoTrans)     trans = 0;
    ಯದಿ (TransA == CblasTrans)       trans = 1;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 0;
    ಯದಿ (TransA == CblasConjTrans)   trans = 1;

    ಯದಿ (Diag == CblasUnit)          unit  = 0;
    ಯದಿ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ಯದಿ (incx == 0)          info =  7;
    ಯದಿ (n < 0)              info =  4;
    ಯದಿ (unit  < 0)          info =  3;
    ಯದಿ (trans < 0)          info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 1;
    ಯದಿ (Uplo == CblasLower)         uplo  = 0;

    ಯದಿ (TransA == CblasNoTrans)     trans = 1;
    ಯದಿ (TransA == CblasTrans)       trans = 0;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 1;
    ಯದಿ (TransA == CblasConjTrans)   trans = 0;

    ಯದಿ (Diag == CblasUnit)          unit  = 0;
    ಯದಿ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ಯದಿ (incx == 0)          info =  7;
    ಯದಿ (n < 0)              info =  4;
    ಯದಿ (unit  < 0)          info =  3;
    ಯದಿ (trans < 0)          info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  ಯದಿ (n == 0) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
