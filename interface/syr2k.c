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

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#ಅನ್ಯಥಾ
#define ERROR_NAME "SSYR2K"
#endif
#ಅನ್ಯಥಾ
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#ಅನ್ಯಥಾ
#define ERROR_NAME "CSYR2K"
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#ಅನ್ಯಥಾ
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, ವರ್ಣ *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ವರ್ಣ uplo_arg  = *UPLO;
  ವರ್ಣ trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  ಪೂರ್ಣಾಙ್ಕ trans;
  ಪೂರ್ಣಾಙ್ಕ nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (ಶೂನ್ಯ *)a;
  args.b = (ಶೂನ್ಯ *)b;
  args.c = (ಶೂನ್ಯ *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
  ಯದಿ (trans_arg == 'C') trans = 1;
#ಅನ್ಯಥಾ
#ifdef HEMM
  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'C') trans = 1;
#ಅನ್ಯಥಾ
  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  ಯದಿ (trans & 1) nrowa = args.k;

  info = 0;

  ಯದಿ (args.ldc < MAX(1,args.n)) info = 12;
  ಯದಿ (args.ldb < MAX(1,nrowa))  info =  9;
  ಯದಿ (args.lda < MAX(1,nrowa))  info =  7;
  ಯದಿ (args.k < 0)               info =  4;
  ಯದಿ (args.n < 0)               info =  3;
  ಯದಿ (trans < 0)                info =  2;
  ಯದಿ (uplo  < 0)                info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_UPLO Uplo, ಗಣನಾ CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint lda,
	   ಶೂನ್ಯ *vb, blasint ldb,
#endif
#ಯದಿ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#ಯದಿ !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ಪೂರ್ಣಾಙ್ಕ uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (ಶೂನ್ಯ *)a;
  args.b = (ಶೂನ್ಯ *)b;
  args.c = (ಶೂನ್ಯ *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (ಶೂನ್ಯ *)&alpha;
#ಅನ್ಯಥಾ
  args.alpha = (ಶೂನ್ಯ *)alpha;
#endif

#ಯದಿ !defined(COMPLEX) || defined(HEMM)
  args.beta  = (ಶೂನ್ಯ *)&beta;
#ಅನ್ಯಥಾ
  args.beta  = (ಶೂನ್ಯ *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper) uplo  = 0;
    ಯದಿ (Uplo == CblasLower) uplo  = 1;

    ಯದಿ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    ಯದಿ (Trans == CblasTrans)       trans = 1;
    ಯದಿ (Trans == CblasConjNoTrans) trans = 0;
    ಯದಿ (Trans == CblasConjTrans)   trans = 1;
#ನೋಚೇದ್ಯದಿ !defined(HEMM)
    ಯದಿ (Trans == CblasTrans)       trans = 1;
#ಅನ್ಯಥಾ
    ಯದಿ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    ಯದಿ (trans & 1) nrowa = args.k;

    ಯದಿ (args.ldc < MAX(1,args.n)) info = 12;
    ಯದಿ (args.ldb < MAX(1,nrowa))  info =  9;
    ಯದಿ (args.lda < MAX(1,nrowa))  info =  7;
    ಯದಿ (args.k < 0)               info =  4;
    ಯದಿ (args.n < 0)               info =  3;
    ಯದಿ (trans < 0)                info =  2;
    ಯದಿ (uplo  < 0)                info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (ಶೂನ್ಯ *)CAlpha;
#endif

    ಯದಿ (Uplo == CblasUpper) uplo  = 1;
    ಯದಿ (Uplo == CblasLower) uplo  = 0;

    ಯದಿ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    ಯದಿ (Trans == CblasTrans)       trans = 0;
    ಯದಿ (Trans == CblasConjNoTrans) trans = 1;
    ಯದಿ (Trans == CblasConjTrans)   trans = 0;
#ನೋಚೇದ್ಯದಿ !defined(HEMM)
    ಯದಿ (Trans == CblasTrans)       trans = 0;
#ಅನ್ಯಥಾ
    ಯದಿ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    ಯದಿ (trans & 1) nrowa = args.k;

    ಯದಿ (args.ldc < MAX(1,args.n)) info = 12;
    ಯದಿ (args.ldb < MAX(1,nrowa))  info =  9;
    ಯದಿ (args.lda < MAX(1,nrowa))  info =  7;
    ಯದಿ (args.k < 0)               info =  4;
    ಯದಿ (args.n < 0)               info =  3;
    ಯದಿ (trans < 0)                info =  2;
    ಯದಿ (uplo  < 0)                info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  ಯದಿ (args.n == 0) ಪ್ರತ್ಯಾಗಮ;

#ifdef DYNAMIC_ARCH
ಬಾಹ್ಯ ವರ್ಣ* gotoblas_corename(ಶೂನ್ಯ);
#endif

#ಯದಿ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಯದಿ defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ಯದಿ defined(DYNAMIC_ARCH)
ಯದಿ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಯದಿ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ಯದಿ (order == CblasRowMajor && n == ldc) {
     ಯದಿ (Trans == CblasNoTrans && k == lda && k == ldb) {
      ಯದಿ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
      }ಅನ್ಯಥಾ ಯದಿ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
      }
     }
     ಅನ್ಯಥಾ ಯದಿ (Trans == CblasTrans && n == lda && n ==ldb) {
      ಯದಿ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
      }ಅನ್ಯಥಾ ಯದಿ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  ಯದಿ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ಅನ್ಯಥಾ {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = ರಿಕ್ತ;
  ಯದಿ (args.n*args.k <1000)
	  args.nthreads =1 ;
  ಅನ್ಯಥಾ
  args.nthreads = ಅಂಕೆ_cpu_avail(3);

  ಯದಿ (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ifdef SMP

  } ಅನ್ಯಥಾ {

    syrk_thread(mode, &args, ರಿಕ್ತ, ರಿಕ್ತ, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
