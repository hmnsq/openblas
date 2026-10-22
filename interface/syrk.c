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

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 109944.
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QSYRK "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DSYRK "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SSYRK "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 14824.
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XSYRK "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZSYRK "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CSYRK "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XHERK "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZHERK "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CHERK "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM_MULTITHREAD_ಮಿತಿ
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_MULTITHREAD_ಮಿತಿ 4
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
};

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *UPLO, ಅಕ್ಷರ *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ಅಕ್ಷರ uplo_arg  = *UPLO;
  ಅಕ್ಷರ trans_arg = *TRANS;

  blas_arg_t ವಾದಗಳು;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ದ್ವಿಗುಣ NNK;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ USE_SIMPLE_THREADED_LEVEL3
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

  blasint info;
  ಪೂರ್ಣಾಂಕ uplo;
  ಪೂರ್ಣಾಂಕ trans;
  ಪೂರ್ಣಾಂಕ nrowa;

  PRINT_DEBUG_NAME;

  ವಾದಗಳು.n = *N;
  ವಾದಗಳು.k = *K;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

  ವಾದಗಳು.lda = *ldA;
  ವಾದಗಳು.ldc = *ldC;

  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ಆದರೆ (uplo_arg  == 'U') uplo  = 0;
  ಆದರೆ (uplo_arg  == 'L') uplo  = 1;


#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ಆದರೆ (trans_arg == 'N') trans = 0;
  ಆದರೆ (trans_arg == 'T') trans = 1;
  ಆದರೆ (trans_arg == 'C') trans = 1;
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ HEMM
  ಆದರೆ (trans_arg == 'N') trans = 0;
  ಆದರೆ (trans_arg == 'C') trans = 1;
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ (trans_arg == 'N') trans = 0;
  ಆದರೆ (trans_arg == 'T') trans = 1;
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ_ಅಂತ್ಯ

  nrowa = ವಾದಗಳು.n;
  ಆದರೆ (trans & 1) nrowa = ವಾದಗಳು.k;

  info = 0;

  ಆದರೆ (ವಾದಗಳು.ldc < MAX(1,ವಾದಗಳು.n)) info = 10;
  ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  7;
  ಆದರೆ (ವಾದಗಳು.k < 0)               info =  4;
  ಆದರೆ (ವಾದಗಳು.n < 0)               info =  3;
  ಆದರೆ (trans < 0)                info =  2;
  ಆದರೆ (uplo  < 0)                info =  1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_UPLO Uplo, ಎಣಿಕೆ CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ಆದರೆ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *valpha,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *va, blasint lda,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *vbeta,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *vc, blasint ldc) {
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ COMPLEX
#ಆದರೆ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#ಆದರೆ_ಅಂತ್ಯ
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#ಆದರೆ_ಅಂತ್ಯ

  blas_arg_t ವಾದಗಳು;
  ಪೂರ್ಣಾಂಕ uplo, trans;
  blasint info, nrowa;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ದ್ವಿಗುಣ NNK;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ USE_SIMPLE_THREADED_LEVEL3
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  ವಾದಗಳು.n = n;
  ವಾದಗಳು.k = k;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

  ವಾದಗಳು.lda = lda;
  ವಾದಗಳು.ldc = ldc;

#ಆದರೆ !defined(COMPLEX) || defined(HEMM)
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)&alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)&beta;
#ಇಲ್ಲದಿದ್ದರೆ
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;
#ಆದರೆ_ಅಂತ್ಯ

  trans = -1;
  uplo  = -1;
  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo  = 0;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo  = 1;

    ಆದರೆ (Trans == CblasNoTrans)     trans = 0;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (Trans == CblasTrans)       trans = 1;
    ಆದರೆ (Trans == CblasConjNoTrans) trans = 0;
    ಆದರೆ (Trans == CblasConjTrans)   trans = 1;
#ಅಥವಾಆದರೆ !defined(HEMM)
    ಆದರೆ (Trans == CblasTrans)       trans = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (Trans == CblasConjTrans)   trans = 1;
#ಆದರೆ_ಅಂತ್ಯ

    info = -1;

    nrowa = ವಾದಗಳು.n;
    ಆದರೆ (trans & 1) nrowa = ವಾದಗಳು.k;

    ಆದರೆ (ವಾದಗಳು.ldc < MAX(1,ವಾದಗಳು.n)) info = 10;
    ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  7;
    ಆದರೆ (ವಾದಗಳು.k < 0)               info =  4;
    ಆದರೆ (ವಾದಗಳು.n < 0)               info =  3;
    ಆದರೆ (trans < 0)                info =  2;
    ಆದರೆ (uplo  < 0)                info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo  = 1;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo  = 0;

    ಆದರೆ (Trans == CblasNoTrans)     trans = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (Trans == CblasTrans)       trans = 0;
    ಆದರೆ (Trans == CblasConjNoTrans) trans = 1;
    ಆದರೆ (Trans == CblasConjTrans)   trans = 0;
#ಅಥವಾಆದರೆ !defined(HEMM)
    ಆದರೆ (Trans == CblasTrans)       trans = 0;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (Trans == CblasConjTrans)   trans = 0;
#ಆದರೆ_ಅಂತ್ಯ

    info = -1;

    nrowa = ವಾದಗಳು.n;
    ಆದರೆ (trans & 1) nrowa = ವಾದಗಳು.k;

    ಆದರೆ (ವಾದಗಳು.ldc < MAX(1,ವಾದಗಳು.n)) info = 10;
    ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  7;
    ಆದರೆ (ವಾದಗಳು.k < 0)               info =  4;
    ಆದರೆ (ವಾದಗಳು.n < 0)               info =  3;
    ಆದರೆ (trans < 0)                info =  2;
    ಆದರೆ (uplo  < 0)                info =  1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

  ಆದರೆ (ವಾದಗಳು.n == 0) ಹಿಂತಿರುಗಿ;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DYNAMIC_ARCH
ಬಾಹ್ಯ ಅಕ್ಷರ* gotoblas_corename(ಶೂನ್ಯ);
#ಆದರೆ_ಅಂತ್ಯ
	
#ಆದರೆ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಆದರೆ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
ಆದರೆ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಆದರೆ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ಆದರೆ_ಅಂತ್ಯ
)
#ಆದರೆ_ಅಂತ್ಯ
  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && n == ldc) {
    ಆದರೆ (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_ಊರ್ಧ್ವ ? SSYRK_ನೇರ_ALPHA_BETA_UN : SSYRK_ನೇರ_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      ಹಿಂತಿರುಗಿ;
    } ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_ಊರ್ಧ್ವ ? SSYRK_ನೇರ_ALPHA_BETA_UT : SSYRK_ನೇರ_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      ಹಿಂತಿರುಗಿ;
    }
  }
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ_ಅಂತ್ಯ


  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(0);
  ಆದರೆ (!ಬಫರ್) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

  sa = (FLOAT *)((BLASLONG)ಬಫರ್ + GEMM_ಸ್ಥಾನಾಂತರ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ಗಾತ್ರ_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ಸ್ಥಾನಾಂತರ_B);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ USE_SIMPLE_THREADED_LEVEL3
  ಆದರೆ (!trans){
    ಬಹುಲಕ |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ಇಲ್ಲದಿದ್ದರೆ {
    ಬಹುಲಕ |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  ಬಹುಲಕ |= (uplo  << BLAS_UPLO_ಸರಿಸು);
#ಆದರೆ_ಅಂತ್ಯ

  ವಾದಗಳು.common = ಖಾಲಿ;

  NNK = (ದ್ವಿಗುಣ)(ವಾದಗಳು.n+1)*(ದ್ವಿಗುಣ)ವಾದಗಳು.n*(ದ್ವಿಗುಣ)ವಾದಗಳು.k;
  ಆದರೆ (NNK <= (SMP_ಮಿತಿ_MIN * GEMM_MULTITHREAD_ಮಿತಿ)) {
  ವಾದಗಳು.nthreads = 1;
  } ಇಲ್ಲದಿದ್ದರೆ {
  ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(3);
  }

  ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (syrk[(uplo << 1) | trans ])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ಇಲ್ಲದಿದ್ದರೆ

    syrk_thread(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, syrk[(uplo << 1) | trans ], sa, sb, ವಾದಗಳು.nthreads);

#ಆದರೆ_ಅಂತ್ಯ

  }
#ಆದರೆ_ಅಂತ್ಯ

 blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE, ವಾದಗಳು.n * ವಾದಗಳು.k + ವಾದಗಳು.n * ವಾದಗಳು.n / 2, ವಾದಗಳು.n * ವಾದಗಳು.n * ವಾದಗಳು.k);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
