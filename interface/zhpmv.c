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

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XHPMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZHPMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CHPMV "
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ  ಪೂರ್ಣಾಂಕ (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#ಇಲ್ಲದಿದ್ದರೆ
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#ಆದರೆ_ಅಂತ್ಯ
};

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ  ಪೂರ್ಣಾಂಕ (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#ಇಲ್ಲದಿದ್ದರೆ
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#ಆದರೆ_ಅಂತ್ಯ
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ಅಕ್ಷರ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  ಪೂರ್ಣಾಂಕ uplo;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ಆದರೆ (uplo_arg  == 'U') uplo  = 0;
  ಆದರೆ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ಆದರೆ (incy == 0)          info =  9;
  ಆದರೆ (incx == 0)          info =  6;
  ಆದರೆ (n < 0)              info =  2;
  ಆದರೆ (uplo  < 0)          info =  1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   ಎಣಿಕೆ CBLAS_UPLO Uplo,
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
  FLOAT *ಬಫರ್;
  ಪೂರ್ಣಾಂಕ uplo;
  blasint info;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)         uplo  = 0;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)         uplo  = 1;

    info = -1;

    ಆದರೆ (incy == 0)          info =  9;
    ಆದರೆ (incx == 0)          info =  6;
    ಆದರೆ (n < 0)              info =  2;
    ಆದರೆ (uplo  < 0)          info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)         uplo  = 3;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)         uplo  = 2;

    info = -1;

    ಆದರೆ (incy == 0)          info =  9;
    ಆದರೆ (incx == 0)          info =  6;
    ಆದರೆ (n < 0)              info =  2;
    ಆದರೆ (uplo  < 0)          info =  1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (n == 0) ಹಿಂತಿರುಗಿ;

  ಆದರೆ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ಖಾಲಿ, 0, ಖಾಲಿ, 0);

  ಆದರೆ ((alpha_r == ZERO) && (alpha_i == ZERO)) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0 ) x -= (n - 1) * incx * 2;
  ಆದರೆ (incy < 0 ) y -= (n - 1) * incy * 2;

  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(1);
  ಆದರೆ (!ಬಫರ್) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  nthreads = ಅಂಕೆ_cpu_avail(2);

  ಆದರೆ (n < 50) nthreads = 1;
  ಆದರೆ (nthreads > 2 && n < 500) nthreads = 2;

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, ಬಫರ್, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(4, n * n / 2 + n, n * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
