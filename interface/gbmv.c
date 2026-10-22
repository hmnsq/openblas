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

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGBMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGBMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGBMV "
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಶೂನ್ಯ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
      qgbmv_n, qgbmv_t,
#ಅಥವಾಆದರೆ defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#ಇಲ್ಲದಿದ್ದರೆ
      sgbmv_n, sgbmv_t,
#ಆದರೆ_ಅಂತ್ಯ
};

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#ಅಥವಾಆದರೆ defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#ಇಲ್ಲದಿದ್ದರೆ
      sgbmv_thread_n, sgbmv_thread_t,
#ಆದರೆ_ಅಂತ್ಯ
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  ಅಕ್ಷರ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ಆದರೆ (trans == 'N') i = 0;
  ಆದರೆ (trans == 'T') i = 1;
  ಆದರೆ (trans == 'R') i = 0;
  ಆದರೆ (trans == 'C') i = 1;

  ಆದರೆ (incy == 0)	 info = 13;
  ಆದರೆ (incx == 0)	 info = 10;
  ಆದರೆ (lda < kl + ku + 1) info = 8;
  ಆದರೆ (kl < 0)		 info = 5;
  ಆದರೆ (ku < 0)		 info = 4;
  ಆದರೆ (n < 0)		 info = 3;
  ಆದರೆ (m < 0)		 info = 2;
  ಆದರೆ (i < 0)		 info = 1;

  trans = i;

  ಆದರೆ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   ಎಣಿಕೆ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *ಬಫರ್;
  blasint lenx, leny, info, t;
  ಪೂರ್ಣಾಂಕ trans;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (TransA == CblasNoTrans)     trans = 0;
    ಆದರೆ (TransA == CblasTrans)       trans = 1;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 0;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ಆದರೆ (incy == 0)	 info = 13;
    ಆದರೆ (incx == 0)	 info = 10;
    ಆದರೆ (lda < kl + ku + 1) info = 8;
    ಆದರೆ (kl < 0)		 info = 5;
    ಆದರೆ (ku < 0)		 info = 4;
    ಆದರೆ (n < 0)		 info = 3;
    ಆದರೆ (m < 0)		 info = 2;
    ಆದರೆ (trans < 0)	 info = 1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (TransA == CblasNoTrans)     trans = 1;
    ಆದರೆ (TransA == CblasTrans)       trans = 0;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 1;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ಆದರೆ (incy == 0)	 info = 13;
    ಆದರೆ (incx == 0)	 info = 10;
    ಆದರೆ (lda < kl + ku + 1) info = 8;
    ಆದರೆ (kl < 0)		 info = 5;
    ಆದರೆ (ku < 0)		 info = 4;
    ಆದರೆ (n < 0)		 info = 3;
    ಆದರೆ (m < 0)		 info = 2;
    ಆದರೆ (trans < 0)	 info = 1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ ((m==0) || (n==0)) ಹಿಂತಿರುಗಿ;

  lenx = n;
  leny = m;
  ಆದರೆ (trans) lenx = m;
  ಆದರೆ (trans) leny = n;

  ಆದರೆ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ಖಾಲಿ, 0, ಖಾಲಿ, 0);

  ಆದರೆ (alpha == ZERO) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0) x -= (lenx-1)*incx;
  ಆದರೆ (incy < 0) y -= (leny-1)*incy;

  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(1);
  ಆದರೆ (!ಬಫರ್) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಆದರೆ (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  ಇಲ್ಲದಿದ್ದರೆ
    nthreads = ಅಂಕೆ_cpu_avail(2);

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

  (gbmv[(ಪೂರ್ಣಾಂಕ)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

    (gbmv_thread[(ಪೂರ್ಣಾಂಕ)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, ಬಫರ್, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(1, m * n / 2 + n, m * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
