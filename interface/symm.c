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
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 65536.
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QSYMM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DSYMM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SSYMM "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 8192.
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM3M
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XSYMM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZSYMM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CSYMM "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XHEMM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZHEMM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CHEMM "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XSYMM3M "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZSYMM3M "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CSYMM3M "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XHEMM3M "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZHEMM3M "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CHEMM3M "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM_MULTITHREAD_ಮಿತಿ
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_MULTITHREAD_ಮಿತಿ 4
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DYNAMIC_ARCH
ಬಾಹ್ಯ ಅಕ್ಷರ* gotoblas_corename(ಶೂನ್ಯ);                                
#ಆದರೆ_ಅಂತ್ಯ                                                               

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ಬಹುಲಕ_2	(BLAS_XDOUBLE | BLAS_REAL)
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ಬಹುಲಕ_2	(BLAS_DOUBLE  | BLAS_REAL)
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ಬಹುಲಕ_2	(BLAS_SINGLE  | BLAS_REAL)
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ಬಹುಲಕ_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ಬಹುಲಕ_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ಬಹುಲಕ_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM3M
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
};

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *SIDE, ಅಕ್ಷರ *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ಅಕ್ಷರ side_arg  = *SIDE;
  ಅಕ್ಷರ uplo_arg  = *UPLO;

  blas_arg_t ವಾದಗಳು;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ಆದರೆ defined(SMP) && !defined(NO_AFFINITY)
  ಪೂರ್ಣಾಂಕ nodes;
#ಆದರೆ_ಅಂತ್ಯ
# ಆದರೆ defined(SMP)
  ದ್ವಿಗುಣ MN;
#ಆದರೆ_ಅಂತ್ಯ
  blasint info;
  ಪೂರ್ಣಾಂಕ side;
  ಪೂರ್ಣಾಂಕ uplo;

  PRINT_DEBUG_NAME;

  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  ಆದರೆ (side_arg  == 'L') side  = 0;
  ಆದರೆ (side_arg  == 'R') side  = 1;

  ಆದರೆ (uplo_arg  == 'U') uplo  = 0;
  ಆದರೆ (uplo_arg  == 'L') uplo  = 1;

  ವಾದಗಳು.m = *M;
  ವಾದಗಳು.n = *N;

  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;
  ವಾದಗಳು.ldc = *ldC;

  info = 0;

  ಆದರೆ (ವಾದಗಳು.ldc < MAX(1, ವಾದಗಳು.m)) info = 12;

  ಆದರೆ (!side) {
    ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
    ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;

    ವಾದಗಳು.lda = *ldA;
    ವಾದಗಳು.ldb = *ldB;

    ಆದರೆ (ವಾದಗಳು.ldb < MAX(1, ವಾದಗಳು.m)) info =  9;
    ಆದರೆ (ವಾದಗಳು.lda < MAX(1, ವಾದಗಳು.m)) info =  7;

  } ಇಲ್ಲದಿದ್ದರೆ {
    ವಾದಗಳು.a = (ಶೂನ್ಯ *)b;
    ವಾದಗಳು.b = (ಶೂನ್ಯ *)a;

    ವಾದಗಳು.lda = *ldB;
    ವಾದಗಳು.ldb = *ldA;

  ಆದರೆ (ವಾದಗಳು.lda < MAX(1, ವಾದಗಳು.m)) info =  9;
  ಆದರೆ (ವಾದಗಳು.ldb < MAX(1, ವಾದಗಳು.n)) info =  7;
  }

  ಆದರೆ (ವಾದಗಳು.n   < 0)              info =  4;
  ಆದರೆ (ವಾದಗಳು.m   < 0)              info =  3;
  ಆದರೆ (uplo     < 0)              info =  2;
  ಆದರೆ (side     < 0)              info =  1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_SIDE Side, ಎಣಿಕೆ CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint lda,
	   ಶೂನ್ಯ *vb, blasint ldb,
	   ಶೂನ್ಯ *vbeta,
	   ಶೂನ್ಯ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#ಆದರೆ_ಅಂತ್ಯ

  blas_arg_t ವಾದಗಳು;
  ಪೂರ್ಣಾಂಕ side, uplo;
  blasint info;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ಆದರೆ defined(SMP) && !defined(NO_AFFINITY)
  ಪೂರ್ಣಾಂಕ nodes;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ defined(SMP)
  ದ್ವಿಗುಣ MN;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)&alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)&beta;
#ಇಲ್ಲದಿದ್ದರೆ
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;
#ಆದರೆ_ಅಂತ್ಯ

  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;
  ವಾದಗಳು.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (Side == Cblas_ಎಡ)  side = 0;
    ಆದರೆ (Side == Cblas_ಬಲ) side = 1;

    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo  = 0;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo  = 1;

    info = -1;

    ವಾದಗಳು.m = m;
    ವಾದಗಳು.n = n;

    ಆದರೆ (ವಾದಗಳು.ldc < MAX(1, ವಾದಗಳು.m)) info = 12;

    ಆದರೆ (!side) {
      ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
      ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;

      ವಾದಗಳು.lda = lda;
      ವಾದಗಳು.ldb = ldb;

      ಆದರೆ (ವಾದಗಳು.ldb < MAX(1, ವಾದಗಳು.m)) info =  9;
      ಆದರೆ (ವಾದಗಳು.lda < MAX(1, ವಾದಗಳು.m)) info =  7;

    } ಇಲ್ಲದಿದ್ದರೆ {
      ವಾದಗಳು.a = (ಶೂನ್ಯ *)b;
      ವಾದಗಳು.b = (ಶೂನ್ಯ *)a;

      ವಾದಗಳು.lda = ldb;
      ವಾದಗಳು.ldb = lda;

      ಆದರೆ (ವಾದಗಳು.lda < MAX(1, ವಾದಗಳು.m)) info =  9;
      ಆದರೆ (ವಾದಗಳು.ldb < MAX(1, ವಾದಗಳು.n)) info =  7;
    }

    ಆದರೆ (ವಾದಗಳು.n   < 0)              info =  4;
    ಆದರೆ (ವಾದಗಳು.m   < 0)              info =  3;
    ಆದರೆ (uplo     < 0)              info =  2;
    ಆದರೆ (side     < 0)              info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (Side == Cblas_ಎಡ)  side = 1;
    ಆದರೆ (Side == Cblas_ಬಲ) side = 0;

    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo  = 1;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo  = 0;

    info = -1;

    ವಾದಗಳು.m = n;
    ವಾದಗಳು.n = m;

    ಆದರೆ (ವಾದಗಳು.ldc < MAX(1, ವಾದಗಳು.m)) info = 12;

    ಆದರೆ (!side) {
      ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
      ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;

      ವಾದಗಳು.lda = lda;
      ವಾದಗಳು.ldb = ldb;

      ಆದರೆ (ವಾದಗಳು.ldb < MAX(1, ವಾದಗಳು.m)) info =  9;
      ಆದರೆ (ವಾದಗಳು.lda < MAX(1, ವಾದಗಳು.m)) info =  7;

    } ಇಲ್ಲದಿದ್ದರೆ {
      ವಾದಗಳು.a = (ಶೂನ್ಯ *)b;
      ವಾದಗಳು.b = (ಶೂನ್ಯ *)a;

      ವಾದಗಳು.lda = ldb;
      ವಾದಗಳು.ldb = lda;

      ಆದರೆ (ವಾದಗಳು.lda < MAX(1, ವಾದಗಳು.m)) info =  9;
      ಆದರೆ (ವಾದಗಳು.ldb < MAX(1, ವಾದಗಳು.n)) info =  7;
    }

    ಆದರೆ (ವಾದಗಳು.n   < 0)              info =  4;
    ಆದರೆ (ವಾದಗಳು.m   < 0)              info =  3;
    ಆದರೆ (uplo     < 0)              info =  2;
    ಆದರೆ (side     < 0)              info =  1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

   ಆದರೆ (ವಾದಗಳು.m == 0 || ವಾದಗಳು.n == 0) ಹಿಂತಿರುಗಿ;
#ಆದರೆ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಆದರೆ defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
ಆದರೆ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಆದರೆ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ಆದರೆ_ಅಂತ್ಯ
)
#ಆದರೆ_ಅಂತ್ಯ
   ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && m == lda && n == ldb && n == ldc)
   {
     ಆದರೆ (Side == Cblas_ಎಡ && Uplo == Cblas_ಊರ್ಧ್ವ) {
       SSYMM_ನೇರ_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ಹಿಂತಿರುಗಿ;
     }
     ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (Side == Cblas_ಎಡ && Uplo == Cblas_ಅಧೋ) {
       SSYMM_ನೇರ_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ಹಿಂತಿರುಗಿ;
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
  ವಾದಗಳು.common = ಖಾಲಿ;
  MN = 2.* (ದ್ವಿಗುಣ) ವಾದಗಳು.m * (ದ್ವಿಗುಣ)ವಾದಗಳು.m * (ದ್ವಿಗುಣ) ವಾದಗಳು.n;
  ಆದರೆ (MN <= (SMP_ಮಿತಿ_MIN * (ದ್ವಿಗುಣ) GEMM_MULTITHREAD_ಮಿತಿ) ) {
	  ವಾದಗಳು.nthreads = 1;
  } ಇಲ್ಲದಿದ್ದರೆ {
	 ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(3);
  }
  ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (symm[(side << 1) | uplo ])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ NO_AFFINITY
    nodes = get_ಅಂಕೆ_nodes();

    ಆದರೆ (nodes > 1) {

      ವಾದಗಳು.nthreads /= nodes;

      gemm_thread_mn(ಬಹುಲಕ_2, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } ಇಲ್ಲದಿದ್ದರೆ {
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ಇಲ್ಲದಿದ್ದರೆ

      GEMM_THREAD(ಬಹುಲಕ_2, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, symm[(side << 1) | uplo ], sa, sb, ವಾದಗಳು.nthreads);

#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ NO_AFFINITY
    }
#ಆದರೆ_ಅಂತ್ಯ

  }
#ಆದರೆ_ಅಂತ್ಯ

 blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE,
		       (!side)? ವಾದಗಳು.m * (ವಾದಗಳು.m / 2 + ವಾದಗಳು.n) : ವಾದಗಳು.n * (ವಾದಗಳು.m + ವಾದಗಳು.n / 2),
		       (!side)? 2 * ವಾದಗಳು.m * ವಾದಗಳು.m * ವಾದಗಳು.n : 2 * ವಾದಗಳು.m * ವಾದಗಳು.n * ವಾದಗಳು.n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
