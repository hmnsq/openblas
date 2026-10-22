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
#define ERROR_NAME "XHPR  "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHPR  "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHPR  "
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*hpr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xhpr_U, xhpr_L, xhpr_V, xhpr_M,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zhpr_U, zhpr_L, zhpr_V, zhpr_M,
#ಅನ್ಯಥಾ
  chpr_U, chpr_L, chpr_V, chpr_M,
#endif
};

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*hpr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  xhpr_thread_U, xhpr_thread_L, xhpr_thread_V, xhpr_thread_M,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zhpr_thread_U, zhpr_thread_L, zhpr_thread_V, zhpr_thread_M,
#ಅನ್ಯಥಾ
  chpr_thread_U, chpr_thread_L, chpr_thread_V, chpr_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  ವರ್ಣ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  FLOAT *buffer;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ಯದಿ (incx == 0)          info =  5;
  ಯದಿ (n < 0)              info =  2;
  ಯದಿ (uplo  < 0)          info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order,
	   ಗಣನಾ CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   ಶೂನ್ಯ  *vx, blasint incx,
	   ಶೂನ್ಯ  *va) {

  FLOAT* x = (FLOAT*) vx;
  FLOAT* a = (FLOAT*) va;

  FLOAT *buffer;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  blasint info;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 0;
    ಯದಿ (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ಯದಿ (incx == 0)          info =  5;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 3;
    ಯದಿ (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    ಯದಿ (incx == 0)          info =  5;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  ಯದಿ (n == 0) ಪ್ರತ್ಯಾಗಮ;

  ಯದಿ (alpha == ZERO) ಪ್ರತ್ಯಾಗಮ;

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
  ಯದಿ (n < 100)
    nthreads = 1;
  ಅನ್ಯಥಾ
    nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (nthreads == 1) {
#endif

  (hpr[uplo])(n, alpha, x, incx, a, buffer);

#ifdef SMP

  } ಅನ್ಯಥಾ {

    (hpr_thread[uplo])(n, alpha, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
