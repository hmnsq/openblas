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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <ctype.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#ଅନ୍ୟଥା
#define ERROR_NAME "SSYR2K"
#endif
#ଅନ୍ୟଥା
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#ଅନ୍ୟଥା
#define ERROR_NAME "CSYR2K"
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#ଅନ୍ୟଥା
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#ଅନ୍ୟଥା
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, ଅକ୍ଷର *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ଅକ୍ଷର uplo_arg  = *UPLO;
  ଅକ୍ଷର trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (ଶୂନ୍ୟ *)a;
  args.b = (ଶୂନ୍ୟ *)b;
  args.c = (ଶୂନ୍ୟ *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (ଶୂନ୍ୟ *)alpha;
  args.beta  = (ଶୂନ୍ୟ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ଯଦି (uplo_arg  == 'U') uplo  = 0;
  ଯଦି (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'T') trans = 1;
  ଯଦି (trans_arg == 'C') trans = 1;
#ଅନ୍ୟଥା
#ifdef HEMM
  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'C') trans = 1;
#ଅନ୍ୟଥା
  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  ଯଦି (trans & 1) nrowa = args.k;

  info = 0;

  ଯଦି (args.ldc < MAX(1,args.n)) info = 12;
  ଯଦି (args.ldb < MAX(1,nrowa))  info =  9;
  ଯଦି (args.lda < MAX(1,nrowa))  info =  7;
  ଯଦି (args.k < 0)               info =  4;
  ଯଦି (args.n < 0)               info =  3;
  ଯଦି (trans < 0)                info =  2;
  ଯଦି (uplo  < 0)                info =  1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_UPLO Uplo, ଗଣନା CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *valpha,
	   ଶୂନ୍ୟ *va, blasint lda,
	   ଶୂନ୍ୟ *vb, blasint ldb,
#endif
#ଯଦି !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#ଯଦି !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (ଶୂନ୍ୟ *)a;
  args.b = (ଶୂନ୍ୟ *)b;
  args.c = (ଶୂନ୍ୟ *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (ଶୂନ୍ୟ *)&alpha;
#ଅନ୍ୟଥା
  args.alpha = (ଶୂନ୍ୟ *)alpha;
#endif

#ଯଦି !defined(COMPLEX) || defined(HEMM)
  args.beta  = (ଶୂନ୍ୟ *)&beta;
#ଅନ୍ୟଥା
  args.beta  = (ଶୂନ୍ୟ *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (Uplo == CblasUpper) uplo  = 0;
    ଯଦି (Uplo == CblasLower) uplo  = 1;

    ଯଦି (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    ଯଦି (Trans == CblasTrans)       trans = 1;
    ଯଦି (Trans == CblasConjNoTrans) trans = 0;
    ଯଦି (Trans == CblasConjTrans)   trans = 1;
#ନଚେତଯଦି !defined(HEMM)
    ଯଦି (Trans == CblasTrans)       trans = 1;
#ଅନ୍ୟଥା
    ଯଦି (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    ଯଦି (trans & 1) nrowa = args.k;

    ଯଦି (args.ldc < MAX(1,args.n)) info = 12;
    ଯଦି (args.ldb < MAX(1,nrowa))  info =  9;
    ଯଦି (args.lda < MAX(1,nrowa))  info =  7;
    ଯଦି (args.k < 0)               info =  4;
    ଯଦି (args.n < 0)               info =  3;
    ଯଦି (trans < 0)                info =  2;
    ଯଦି (uplo  < 0)                info =  1;
  }

  ଯଦି (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (ଶୂନ୍ୟ *)CAlpha;
#endif

    ଯଦି (Uplo == CblasUpper) uplo  = 1;
    ଯଦି (Uplo == CblasLower) uplo  = 0;

    ଯଦି (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    ଯଦି (Trans == CblasTrans)       trans = 0;
    ଯଦି (Trans == CblasConjNoTrans) trans = 1;
    ଯଦି (Trans == CblasConjTrans)   trans = 0;
#ନଚେତଯଦି !defined(HEMM)
    ଯଦି (Trans == CblasTrans)       trans = 0;
#ଅନ୍ୟଥା
    ଯଦି (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    ଯଦି (trans & 1) nrowa = args.k;

    ଯଦି (args.ldc < MAX(1,args.n)) info = 12;
    ଯଦି (args.ldb < MAX(1,nrowa))  info =  9;
    ଯଦି (args.lda < MAX(1,nrowa))  info =  7;
    ଯଦି (args.k < 0)               info =  4;
    ଯଦି (args.n < 0)               info =  3;
    ଯଦି (trans < 0)                info =  2;
    ଯଦି (uplo  < 0)                info =  1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

  ଯଦି (args.n == 0) ଫେରନ୍ତୁ;

#ifdef DYNAMIC_ARCH
ବାହ୍ୟ ଅକ୍ଷର* gotoblas_corename(ଶୂନ୍ୟ);
#endif

#ଯଦି !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ଯଦି defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ଯଦି defined(DYNAMIC_ARCH)
ଯଦି (strcmp(gotoblas_corename(), "armv9sme") == 0
#ଯଦି defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ଯଦି (order == CblasRowMajor && n == ldc) {
     ଯଦି (Trans == CblasNoTrans && k == lda && k == ldb) {
      ଯଦି (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ଫେରନ୍ତୁ;
      }ଅନ୍ୟଥା ଯଦି (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ଫେରନ୍ତୁ;
      }
     }
     ଅନ୍ୟଥା ଯଦି (Trans == CblasTrans && n == lda && n ==ldb) {
      ଯଦି (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ଫେରନ୍ତୁ;
      }ଅନ୍ୟଥା ଯଦି (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ଫେରନ୍ତୁ;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ଯଦି (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  ଯଦି (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ଅନ୍ୟଥା {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = ଖାଲି;
  ଯଦି (args.n*args.k <1000)
	  args.nthreads =1 ;
  ଅନ୍ୟଥା
  args.nthreads = ଅଙ୍କ_cpu_avail(3);

  ଯଦି (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP

  } ଅନ୍ୟଥା {

    syrk_thread(mode, &args, ଖାଲି, ଖାଲି, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
