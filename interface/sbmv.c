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
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QSBMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DSBMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SSBMV "
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ  ಪೂರ್ಣಾಂಕ (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  qsbmv_U, qsbmv_L,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  dsbmv_U, dsbmv_L,
#ಇಲ್ಲದಿದ್ದರೆ
  ssbmv_U, ssbmv_L,
#ಆದರೆ_ಅಂತ್ಯ
};

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
ಸ್ಥಿರ  ಪೂರ್ಣಾಂಕ (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  qsbmv_thread_U, qsbmv_thread_L,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  dsbmv_thread_U, dsbmv_thread_L,
#ಇಲ್ಲದಿದ್ದರೆ
  ssbmv_thread_U, ssbmv_thread_L,
#ಆದರೆ_ಅಂತ್ಯ
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ಅಕ್ಷರ uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  ಪೂರ್ಣಾಂಕ uplo;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ಆದರೆ (uplo_arg  == 'U') uplo  = 0;
  ಆದರೆ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ಆದರೆ (incy == 0)          info = 11;
  ಆದರೆ (incx == 0)          info =  8;
  ಆದರೆ (lda  < k + 1)       info =  6;
  ಆದರೆ (k < 0)              info =  3;
  ಆದರೆ (n < 0)              info =  2;
  ಆದರೆ (uplo  < 0)          info =  1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   ಎಣಿಕೆ CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *ಬಫರ್;
  ಪೂರ್ಣಾಂಕ uplo;
  blasint info;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)         uplo  = 0;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)         uplo  = 1;

    info = -1;

    ಆದರೆ (incy == 0)          info = 11;
    ಆದರೆ (incx == 0)          info =  8;
    ಆದರೆ (lda  < k + 1)       info =  6;
    ಆದರೆ (k < 0)              info =  3;
    ಆದರೆ (n < 0)              info =  2;
    ಆದರೆ (uplo  < 0)          info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)         uplo  = 1;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)         uplo  = 0;

    info = -1;

    ಆದರೆ (incy == 0)          info = 11;
    ಆದರೆ (incx == 0)          info =  8;
    ಆದರೆ (lda  < k + 1)       info =  6;
    ಆದರೆ (k < 0)              info =  3;
    ಆದರೆ (n < 0)              info =  2;
    ಆದರೆ (uplo  < 0)          info =  1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (n == 0) ಹಿಂತಿರುಗಿ;

  ಆದರೆ (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), ಖಾಲಿ, 0, ಖಾಲಿ, 0);

  ಆದರೆ (alpha == ZERO) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0 ) x -= (n - 1) * incx;
  ಆದರೆ (incy < 0 ) y -= (n - 1) * incy;

  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(1);
  ಆದರೆ (!ಬಫರ್) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

  (sbmv[uplo])(n, k, alpha, a, lda, x, incx, y, incy, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  } ಇಲ್ಲದಿದ್ದರೆ {

    (sbmv_thread[uplo])(n, k, alpha, a, lda, x, incx, y, incy, ಬಫರ್, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(1, n * k / 2 + n, n * k);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
