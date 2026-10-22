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

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ifdef XDOUBLE
#define ERROR_NAME "XSBMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZSBMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CSBMV "
#endif

ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
  xsbmv_U, xsbmv_L,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zsbmv_U, zsbmv_L,
#ಅನ್ಯಥಾ
  csbmv_U, csbmv_L,
#endif
};

#ifdef SMPTEST
ಸ್ಥಿರ  ಪೂರ್ಣಾಙ್ಕ (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  xsbmv_thread_U, xsbmv_thread_L,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  zsbmv_thread_U, zsbmv_thread_L,
#ಅನ್ಯಥಾ
  csbmv_thread_U, csbmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT
ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *b, blasint *INCX, FLOAT *BETA, FLOAT *c, blasint *INCY){

  ವರ್ಣ uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
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

  ಯದಿ (incy == 0)          info = 11;
  ಯದಿ (incx == 0)          info =  8;
  ಯದಿ (lda  < k + 1)       info =  6;
  ಯದಿ (k < 0)              info =  3;
  ಯದಿ (n < 0)              info =  2;
  ಯದಿ (uplo  < 0)          info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  ಯದಿ (n == 0) ಪ್ರತ್ಯಾಗಮ;

  ಯದಿ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, c, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ ((alpha_r == ZERO) && (alpha_i == ZERO)) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0 ) b -= (n - 1) * incx * COMPSIZE;
  ಯದಿ (incy < 0 ) c -= (n - 1) * incy * COMPSIZE;

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

  (sbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, b, incx, c, incy, buffer);

#ifdef SMPTEST
  } ಅನ್ಯಥಾ {

    (sbmv_thread[uplo])(n, k, ALPHA, a, lda, b, incx, c, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
