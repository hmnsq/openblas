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
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QTBMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DTBMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "STBMV "
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#ಇಲ್ಲದಿದ್ದರೆ
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#ಆದರೆ_ಅಂತ್ಯ
};

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#ಅಥವಾಆದರೆ defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#ಇಲ್ಲದಿದ್ದರೆ
  stbmv_thread_NUU, stbmv_thread_NUN, stbmv_thread_NLU, stbmv_thread_NLN,
  stbmv_thread_TUU, stbmv_thread_TUN, stbmv_thread_TLU, stbmv_thread_TLN,
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
  ಆದರೆ (trans_arg == 'R') trans = 0;
  ಆದರೆ (trans_arg == 'C') trans = 1;

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
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

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
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 0;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 1;

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
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 1;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 0;

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

  ಆದರೆ (incx < 0 ) x -= (n - 1) * incx;

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

  FUNCTION_PROFILE_ಕೊನೆ(1, n * k / 2 + n, n * k);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
