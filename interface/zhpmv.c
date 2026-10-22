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
#define ERROR_NAME "XHPMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHPMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHPMV "
#endif

ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#ಅನ್ಯಥಾ
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#endif
};

#ifdef SMP
ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#ಅನ್ಯಥಾ
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ವರ್ಣ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

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

  ಯದಿ (incy == 0)          info =  9;
  ಯದಿ (incx == 0)          info =  6;
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
	   ಶೂನ್ಯ *VALPHA,
	   ಶೂನ್ಯ  *va,
	   ಶೂನ್ಯ  *vx, blasint incx,
	   ಶೂನ್ಯ *VBETA,
	   ಶೂನ್ಯ  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
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

    ಯದಿ (incy == 0)          info =  9;
    ಯದಿ (incx == 0)          info =  6;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 3;
    ಯದಿ (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    ಯದಿ (incy == 0)          info =  9;
    ಯದಿ (incx == 0)          info =  6;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  ಯದಿ (n == 0) ಪ್ರತ್ಯಾಗಮ;

  ಯದಿ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ ((alpha_r == ZERO) && (alpha_i == ZERO)) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0 ) x -= (n - 1) * incx * 2;
  ಯದಿ (incy < 0 ) y -= (n - 1) * incy * 2;

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

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, buffer);

#ifdef SMP
  } ಅನ್ಯಥಾ {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
