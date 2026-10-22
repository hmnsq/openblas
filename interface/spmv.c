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
#define ERROR_NAME "QSPMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DSPMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SSPMV "
#endif

ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*spmv[])(BLASLONG, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  qspmv_U, qspmv_L,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dspmv_U, dspmv_L,
#ಅನ್ಯಥಾ
  sspmv_U, sspmv_L,
#endif
};

#ifdef SMPTEST
ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*spmv_thread[])(BLASLONG, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  qspmv_thread_U, qspmv_thread_L,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dspmv_thread_U, dspmv_thread_L,
#ಅನ್ಯಥಾ
  sspmv_thread_U, sspmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ವರ್ಣ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  FLOAT *buffer;
#ifdef SMPTEST
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
	   FLOAT alpha,
	   FLOAT  *a,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  blasint info;
#ifdef SMPTEST
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
    ಯದಿ (Uplo == CblasUpper)         uplo  = 1;
    ಯದಿ (Uplo == CblasLower)         uplo  = 0;

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

  ಯದಿ (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ (alpha == ZERO) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0 ) x -= (n - 1) * incx;
  ಯದಿ (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ifdef SMPTEST
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (nthreads == 1) {
#endif

  (spmv[uplo])(n, alpha, a, x, incx, y, incy, buffer);

#ifdef SMPTEST
  } ಅನ್ಯಥಾ {

    (spmv_thread[uplo])(n, alpha, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
