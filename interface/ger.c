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
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ __64BIT__
#ವ್ಯಾಖ್ಯಾನಿಸು SMPTEST 1
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGER  "
#ಅಥವಾಆದರೆ defined DOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGER  "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGER  "
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿಸು GER		GERU_K

#ಆದರೆ   defined XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು GER_THREAD	qger_thread
#ಅಥವಾಆದರೆ defined DOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು GER_THREAD	dger_thread
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು GER_THREAD	sger_thread
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(blasint *M, blasint *N, FLOAT *Alpha,
	  FLOAT *x, blasint *INCX,
	  FLOAT *y, blasint *INCY,
	  FLOAT *a, blasint *LDA){

  blasint    m     = *M;
  blasint    n     = *N;
  FLOAT  alpha = *Alpha;
  blasint    incx  = *INCX;
  blasint    incy  = *INCY;
  blasint    lda   = *LDA;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  blasint info;

  PRINT_DEBUG_NAME;

  info = 0;

  ಆದರೆ (lda < MAX(1,m)) info = 9;
  ಆದರೆ (incy == 0)      info = 7;
  ಆದರೆ (incx == 0)      info = 5;
  ಆದರೆ (n < 0)          info = 2;
  ಆದರೆ (m < 0)          info = 1;

  ಆದರೆ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a, blasint lda) {

  FLOAT *ಬಫರ್;
  blasint info, t;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    info = -1;

    ಆದರೆ (lda < MAX(1,m)) info = 9;
    ಆದರೆ (incy == 0)      info = 7;
    ಆದರೆ (incx == 0)      info = 5;
    ಆದರೆ (n < 0)          info = 2;
    ಆದರೆ (m < 0)          info = 1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

    t    = incx;
    incx = incy;
    incy = t;

    ಬಫರ್ = x;
    x = y;
    y = ಬಫರ್;

    ಆದರೆ (lda < MAX(1,m)) info = 9;
    ಆದರೆ (incy == 0)      info = 7;
    ಆದರೆ (incx == 0)      info = 5;
    ಆದರೆ (n < 0)          info = 2;
    ಆದರೆ (m < 0)          info = 1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ

  /*     Quick return if possible. */
  ಆದರೆ (m == 0 || n == 0) ಹಿಂತಿರುಗಿ;
  ಆದರೆ (alpha == 0.) ಹಿಂತಿರುಗಿ;

  ಆದರೆ (incx == 1 && incy == 1 && 1L*m*n <= 2048 *GEMM_MULTITHREAD_ಮಿತಿ) {
    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, ಖಾಲಿ);
    ಹಿಂತಿರುಗಿ;
  }  

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incy < 0) y -= (n - 1) * incy;
  ಆದರೆ (incx < 0) x -= (m - 1) * incx;

  ರಾಶಿ_ALLOC(m, FLOAT, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  // Threshold chosen so that speed-up is > 1 on a Xeon E5-2630
  ಆದರೆ(1L * m * n > 20480L * GEMM_MULTITHREAD_ಮಿತಿ)
    nthreads = ಅಂಕೆ_cpu_avail(2);
  ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ(1L * m * n > 2048L * GEMM_MULTITHREAD_ಮಿತಿ)
	  nthreads = MIN(3,ಅಂಕೆ_cpu_avail(2));
  ಇಲ್ಲದಿದ್ದರೆ
    nthreads = 1;

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    GER(m, n, 0, alpha, x, incx, y, incy, a, lda, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMPTEST
  } ಇಲ್ಲದಿದ್ದರೆ {

    GER_THREAD(m, n, alpha, x, incx, y, incy, a, lda, ಬಫರ್, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  ರಾಶಿ_FREE(ಬಫರ್);
  FUNCTION_PROFILE_ಕೊನೆ(1, m * n + m + n, 2 * m * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
