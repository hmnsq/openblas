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
#define ERROR_NAME "XTPMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZTPMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CTPMV "
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  xtpmv_NUU, xtpmv_NUN, xtpmv_NLU, xtpmv_NLN,
  xtpmv_TUU, xtpmv_TUN, xtpmv_TLU, xtpmv_TLN,
  xtpmv_RUU, xtpmv_RUN, xtpmv_RLU, xtpmv_RLN,
  xtpmv_CUU, xtpmv_CUN, xtpmv_CLU, xtpmv_CLN,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ztpmv_NUU, ztpmv_NUN, ztpmv_NLU, ztpmv_NLN,
  ztpmv_TUU, ztpmv_TUN, ztpmv_TLU, ztpmv_TLN,
  ztpmv_RUU, ztpmv_RUN, ztpmv_RLU, ztpmv_RLN,
  ztpmv_CUU, ztpmv_CUN, ztpmv_CLU, ztpmv_CLN,
#ಅನ್ಯಥಾ
  ctpmv_NUU, ctpmv_NUN, ctpmv_NLU, ctpmv_NLN,
  ctpmv_TUU, ctpmv_TUN, ctpmv_TLU, ctpmv_TLN,
  ctpmv_RUU, ctpmv_RUN, ctpmv_RLU, ctpmv_RLN,
  ctpmv_CUU, ctpmv_CUN, ctpmv_CLU, ctpmv_CLN,
#endif
};

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  xtpmv_thread_NUU, xtpmv_thread_NUN, xtpmv_thread_NLU, xtpmv_thread_NLN,
  xtpmv_thread_TUU, xtpmv_thread_TUN, xtpmv_thread_TLU, xtpmv_thread_TLN,
  xtpmv_thread_RUU, xtpmv_thread_RUN, xtpmv_thread_RLU, xtpmv_thread_RLN,
  xtpmv_thread_CUU, xtpmv_thread_CUN, xtpmv_thread_CLU, xtpmv_thread_CLN,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ztpmv_thread_NUU, ztpmv_thread_NUN, ztpmv_thread_NLU, ztpmv_thread_NLN,
  ztpmv_thread_TUU, ztpmv_thread_TUN, ztpmv_thread_TLU, ztpmv_thread_TLN,
  ztpmv_thread_RUU, ztpmv_thread_RUN, ztpmv_thread_RLU, ztpmv_thread_RLN,
  ztpmv_thread_CUU, ztpmv_thread_CUN, ztpmv_thread_CLU, ztpmv_thread_CLN,
#ಅನ್ಯಥಾ
  ctpmv_thread_NUU, ctpmv_thread_NUN, ctpmv_thread_NLU, ctpmv_thread_NLN,
  ctpmv_thread_TUU, ctpmv_thread_TUN, ctpmv_thread_TLU, ctpmv_thread_TLN,
  ctpmv_thread_RUU, ctpmv_thread_RUN, ctpmv_thread_RLU, ctpmv_thread_RLN,
  ctpmv_thread_CUU, ctpmv_thread_CUN, ctpmv_thread_CLU, ctpmv_thread_CLN,
#endif
};
#endif

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
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
  ಯದಿ (trans_arg == 'R') trans = 2;
  ಯದಿ (trans_arg == 'C') trans = 3;

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
	   blasint n, ಶೂನ್ಯ  *va, ಶೂನ್ಯ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  ಪೂರ್ಣಾಙ್ಕ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 0;
    ಯದಿ (Uplo == CblasLower)         uplo  = 1;

    ಯದಿ (TransA == CblasNoTrans)     trans = 0;
    ಯದಿ (TransA == CblasTrans)       trans = 1;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 2;
    ಯದಿ (TransA == CblasConjTrans)   trans = 3;

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
    ಯದಿ (TransA == CblasConjNoTrans) trans = 3;
    ಯದಿ (TransA == CblasConjTrans)   trans = 2;

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

  ಯದಿ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ifdef SMP
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (n < 50) nthreads = 1;
  ಯದಿ (nthreads > 2 && n < 500) nthreads = 2;

  ಯದಿ (nthreads == 1) {
#endif

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

#ifdef SMP

  } ಅನ್ಯಥಾ {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
