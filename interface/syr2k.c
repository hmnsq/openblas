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
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QSYR2K"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DSYR2K"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SSYR2K"
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XSYR2K"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZSYR2K"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CSYR2K"
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XHER2K"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZHER2K"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CHER2K"
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#ಇಲ್ಲದಿದ್ದರೆ
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#ಆದರೆ_ಅಂತ್ಯ
};

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *UPLO, ಅಕ್ಷರ *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ಅಕ್ಷರ uplo_arg  = *UPLO;
  ಅಕ್ಷರ trans_arg = *TRANS;

  blas_arg_t ವಾದಗಳು;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
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

  blasint info;
  ಪೂರ್ಣಾಂಕ uplo;
  ಪೂರ್ಣಾಂಕ trans;
  ಪೂರ್ಣಾಂಕ nrowa;

  PRINT_DEBUG_NAME;

  ವಾದಗಳು.n = *N;
  ವಾದಗಳು.k = *K;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;
  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

  ವಾದಗಳು.lda = *ldA;
  ವಾದಗಳು.ldb = *ldB;
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

  ಆದರೆ (ವಾದಗಳು.ldc < MAX(1,ವಾದಗಳು.n)) info = 12;
  ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,nrowa))  info =  9;
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
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint lda,
	   ಶೂನ್ಯ *vb, blasint ldb,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *vbeta,
#ಆದರೆ_ಅಂತ್ಯ
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT *c,
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *vc,
#ಆದರೆ_ಅಂತ್ಯ
           blasint ldc) {

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#ಆದರೆ !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#ಆದರೆ_ಅಂತ್ಯ
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#ಆದರೆ_ಅಂತ್ಯ

  blas_arg_t ವಾದಗಳು;
  ಪೂರ್ಣಾಂಕ uplo, trans;
  blasint info, nrowa;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ HEMM
  FLOAT CAlpha[2];
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
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

  PRINT_DEBUG_CNAME;

  ವಾದಗಳು.n = n;
  ವಾದಗಳು.k = k;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;
  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

  ವಾದಗಳು.lda = lda;
  ವಾದಗಳು.ldb = ldb;
  ವಾದಗಳು.ldc = ldc;

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)&alpha;
#ಇಲ್ಲದಿದ್ದರೆ
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ !defined(COMPLEX) || defined(HEMM)
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)&beta;
#ಇಲ್ಲದಿದ್ದರೆ
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

    ಆದರೆ (ವಾದಗಳು.ldc < MAX(1,ವಾದಗಳು.n)) info = 12;
    ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,nrowa))  info =  9;
    ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  7;
    ಆದರೆ (ವಾದಗಳು.k < 0)               info =  4;
    ಆದರೆ (ವಾದಗಳು.n < 0)               info =  3;
    ಆದರೆ (trans < 0)                info =  2;
    ಆದರೆ (uplo  < 0)                info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)CAlpha;
#ಆದರೆ_ಅಂತ್ಯ

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

    ಆದರೆ (ವಾದಗಳು.ldc < MAX(1,ವಾದಗಳು.n)) info = 12;
    ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,nrowa))  info =  9;
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
#ಆದರೆ defined(ARCH_ARM64) && (defined(USE_SSYR2_K_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
ಆದರೆ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಆದರೆ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ಆದರೆ_ಅಂತ್ಯ
)
#ಆದರೆ_ಅಂತ್ಯ
   ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && n == ldc) {
     ಆದರೆ (Trans == CblasNoTrans && k == lda && k == ldb) {
      ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) {
        SSYR2_K_ನೇರ_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಹಿಂತಿರುಗಿ;
      }ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (Uplo == Cblas_ಅಧೋ) {
        SSYR2_K_ನೇರ_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಹಿಂತಿರುಗಿ;
      }
     }
     ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (Trans == CblasTrans && n == lda && n ==ldb) {
      ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) {
        SSYR2_K_ನೇರ_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಹಿಂತಿರುಗಿ;
      }ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (Uplo == Cblas_ಅಧೋ) {
        SSYR2_K_ನೇರ_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಹಿಂತಿರುಗಿ;
      }
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
  ಆದರೆ (!trans){
    ಬಹುಲಕ |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ಇಲ್ಲದಿದ್ದರೆ {
    ಬಹುಲಕ |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  ಬಹುಲಕ |= (uplo  << BLAS_UPLO_ಸರಿಸು);

  ವಾದಗಳು.common = ಖಾಲಿ;
  ಆದರೆ (ವಾದಗಳು.n*ವಾದಗಳು.k <1000)
	  ವಾದಗಳು.nthreads =1 ;
  ಇಲ್ಲದಿದ್ದರೆ
  ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(3);

  ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (syr2k[(uplo << 1) | trans ])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP

  } ಇಲ್ಲದಿದ್ದರೆ {

    syrk_thread(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, syr2k[(uplo << 1) | trans ], sa, sb, ವಾದಗಳು.nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE, 2 * ವಾದಗಳು.n * ವಾದಗಳು.k + ವಾದಗಳು.n * ವಾದಗಳು.n, 2 * ವಾದಗಳು.n * ವಾದಗಳು.n * ವಾದಗಳು.k);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
