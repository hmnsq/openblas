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
#define ERROR_NAME "QSPR2 "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DSPR2 "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SSPR2 "
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*spr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  qspr2_U, qspr2_L,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dspr2_U, dspr2_L,
#ಅನ್ಯಥಾ
  sspr2_U, sspr2_L,
#endif
};

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*spr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  qspr2_thread_U, qspr2_thread_L,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  dspr2_thread_U, dspr2_thread_L,
#ಅನ್ಯಥಾ
  sspr2_thread_U, sspr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  ವರ್ಣ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
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

  ಯದಿ (incy == 0)          info =  7;
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
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a) {

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

    ಯದಿ (incy == 0)          info =  7;
    ಯದಿ (incx == 0)          info =  5;
    ಯದಿ (n < 0)              info =  2;
    ಯದಿ (uplo  < 0)          info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper)         uplo  = 1;
    ಯದಿ (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    ಯದಿ (incy == 0)          info =  7;
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

  ಯದಿ (incx == 1 && incy == 1 && n < 50) {
    blasint i;
    ಯದಿ (!uplo) {
      ಪ್ರತಿ (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, ರಿಕ್ತ, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, ರಿಕ್ತ, 0);
        a += i + 1;
      }
    } ಅನ್ಯಥಾ {
      ಪ್ರತಿ (i = 0; i < n; i++){
	AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, ರಿಕ್ತ, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, ರಿಕ್ತ, 0);
        a += n - i;
      }
    }
    ಪ್ರತ್ಯಾಗಮ;
  }

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

#ifdef SMP
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (n < 50) nthreads = 1;
  ಯದಿ (nthreads > 2 && n < 500) nthreads = 2;

  ಯದಿ (nthreads == 1) {
#endif

  (spr2[uplo])(n, alpha, x, incx, y, incy, a, buffer);

#ifdef SMP
  } ಅನ್ಯಥಾ {

    (spr2_thread[uplo])(n, alpha, x, incx, y, incy, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
