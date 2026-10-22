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
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XTBMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZTBMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CTBMV "
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  xtbmv_NUU, xtbmv_NUN, xtbmv_NLU, xtbmv_NLN,
  xtbmv_TUU, xtbmv_TUN, xtbmv_TLU, xtbmv_TLN,
  xtbmv_RUU, xtbmv_ಸತತ, xtbmv_RLU, xtbmv_RLN,
  xtbmv_CUU, xtbmv_CUN, xtbmv_CLU, xtbmv_CLN,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ztbmv_NUU, ztbmv_NUN, ztbmv_NLU, ztbmv_NLN,
  ztbmv_TUU, ztbmv_TUN, ztbmv_TLU, ztbmv_TLN,
  ztbmv_RUU, ztbmv_ಸತತ, ztbmv_RLU, ztbmv_RLN,
  ztbmv_CUU, ztbmv_CUN, ztbmv_CLU, ztbmv_CLN,
#ಇಲ್ಲದಿದ್ದರೆ
  ctbmv_NUU, ctbmv_NUN, ctbmv_NLU, ctbmv_NLN,
  ctbmv_TUU, ctbmv_TUN, ctbmv_TLU, ctbmv_TLN,
  ctbmv_RUU, ctbmv_ಸತತ, ctbmv_RLU, ctbmv_RLN,
  ctbmv_CUU, ctbmv_CUN, ctbmv_CLU, ctbmv_CLN,
#ಆದರೆ_ಅಂತ್ಯ
};

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  xtbmv_thread_NUU, xtbmv_thread_NUN, xtbmv_thread_NLU, xtbmv_thread_NLN,
  xtbmv_thread_TUU, xtbmv_thread_TUN, xtbmv_thread_TLU, xtbmv_thread_TLN,
  xtbmv_thread_RUU, xtbmv_thread_ಸತತ, xtbmv_thread_RLU, xtbmv_thread_RLN,
  xtbmv_thread_CUU, xtbmv_thread_CUN, xtbmv_thread_CLU, xtbmv_thread_CLN,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ztbmv_thread_NUU, ztbmv_thread_NUN, ztbmv_thread_NLU, ztbmv_thread_NLN,
  ztbmv_thread_TUU, ztbmv_thread_TUN, ztbmv_thread_TLU, ztbmv_thread_TLN,
  ztbmv_thread_RUU, ztbmv_thread_ಸತತ, ztbmv_thread_RLU, ztbmv_thread_RLN,
  ztbmv_thread_CUU, ztbmv_thread_CUN, ztbmv_thread_CLU, ztbmv_thread_CLN,
#ಇಲ್ಲದಿದ್ದರೆ
  ctbmv_thread_NUU, ctbmv_thread_NUN, ctbmv_thread_NLU, ctbmv_thread_NLN,
  ctbmv_thread_TUU, ctbmv_thread_TUN, ctbmv_thread_TLU, ctbmv_thread_TLN,
  ctbmv_thread_RUU, ctbmv_thread_ಸತತ, ctbmv_thread_RLU, ctbmv_thread_RLN,
  ctbmv_thread_CUU, ctbmv_thread_CUN, ctbmv_thread_CLU, ctbmv_thread_CLN,
#ಆದರೆ_ಅಂತ್ಯ
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *UPLO, ಅಕ್ಷರ *TRANS, ಅಕ್ಷರ *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  ಅಕ್ಷರ uplo_arg  = *UPLO;
  ಅಕ್ಷರ trans_arg = *TRANS;
  ಅಕ್ಷರ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  ಪೂರ್ಣಾಂಕ uplo;
  ಪೂರ್ಣಾಂಕ unit;
  ಪೂರ್ಣಾಂಕ trans;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ಆದರೆ (trans_arg == 'N') trans = 0;
  ಆದರೆ (trans_arg == 'T') trans = 1;
  ಆದರೆ (trans_arg == 'R') trans = 2;
  ಆದರೆ (trans_arg == 'C') trans = 3;

  ಆದರೆ (diag_arg  == 'U') unit  = 0;
  ಆದರೆ (diag_arg  == 'N') unit  = 1;

  ಆದರೆ (uplo_arg  == 'U') uplo  = 0;
  ಆದರೆ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ಆದರೆ (incx == 0)          info =  9;
  ಆದರೆ (lda < k + 1)        info =  7;
  ಆದರೆ (k < 0)              info =  5;
  ಆದರೆ (n < 0)              info =  4;
  ಆದರೆ (unit  < 0)          info =  3;
  ಆದರೆ (trans < 0)          info =  2;
  ಆದರೆ (uplo  < 0)          info =  1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_UPLO Uplo,
	   ಎಣಿಕೆ CBLAS_TRANSPOSE TransA, ಎಣಿಕೆ CBLAS_DIAG Diag,
	   blasint n, blasint k, ಶೂನ್ಯ  *va, blasint lda, ಶೂನ್ಯ  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;

  ಪೂರ್ಣಾಂಕ trans, uplo, unit;
  blasint info;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)         uplo  = 0;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)         uplo  = 1;

    ಆದರೆ (TransA == CblasNoTrans)     trans = 0;
    ಆದರೆ (TransA == CblasTrans)       trans = 1;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 2;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 3;

    ಆದರೆ (Diag == CblasUnit)          unit  = 0;
    ಆದರೆ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ಆದರೆ (incx == 0)          info =  9;
    ಆದರೆ (lda < k + 1)        info =  7;
    ಆದರೆ (k < 0)              info =  5;
    ಆದರೆ (n < 0)              info =  4;
    ಆದರೆ (unit  < 0)          info =  3;
    ಆದರೆ (trans < 0)          info =  2;
    ಆದರೆ (uplo  < 0)          info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)         uplo  = 1;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)         uplo  = 0;

    ಆದರೆ (TransA == CblasNoTrans)     trans = 1;
    ಆದರೆ (TransA == CblasTrans)       trans = 0;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 3;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 2;

    ಆದರೆ (Diag == CblasUnit)          unit  = 0;
    ಆದರೆ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ಆದರೆ (incx == 0)          info =  9;
    ಆದರೆ (lda < k + 1)        info =  7;
    ಆದರೆ (k < 0)              info =  5;
    ಆದರೆ (n < 0)              info =  4;
    ಆದರೆ (unit  < 0)          info =  3;
    ಆದರೆ (trans < 0)          info =  2;
    ಆದರೆ (uplo  < 0)          info =  1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (n == 0) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0 ) x -= (n - 1) * incx * 2;

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

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, ಬಫರ್, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(4, n * k / 2 + n, n * k);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
