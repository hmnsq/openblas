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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <ctype.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#இல்லையெனில்
#define ERROR_NAME "SSYR2K"
#endif
#இல்லையெனில்
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#இல்லையெனில்
#define ERROR_NAME "CSYR2K"
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#இல்லையெனில்
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

நிலையான முழு (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#இல்லையெனில்
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, எழுத்து *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  எழுத்து uplo_arg  = *UPLO;
  எழுத்து trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  முழு uplo;
  முழு trans;
  முழு nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (வெறுமை *)a;
  args.b = (வெறுமை *)b;
  args.c = (வெறுமை *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'T') trans = 1;
  எனில் (trans_arg == 'C') trans = 1;
#இல்லையெனில்
#ifdef HEMM
  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'C') trans = 1;
#இல்லையெனில்
  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  எனில் (trans & 1) nrowa = args.k;

  info = 0;

  எனில் (args.ldc < MAX(1,args.n)) info = 12;
  எனில் (args.ldb < MAX(1,nrowa))  info =  9;
  எனில் (args.lda < MAX(1,nrowa))  info =  7;
  எனில் (args.k < 0)               info =  4;
  எனில் (args.n < 0)               info =  3;
  எனில் (trans < 0)                info =  2;
  எனில் (uplo  < 0)                info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_UPLO Uplo, எண்ணிக்கை CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#இல்லையெனில்
	   வெறுமை *valpha,
	   வெறுமை *va, blasint lda,
	   வெறுமை *vb, blasint ldb,
#endif
#எனில் !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#இல்லையெனில்
	   வெறுமை *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#இல்லையெனில்
	   வெறுமை *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#எனில் !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  முழு uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (வெறுமை *)a;
  args.b = (வெறுமை *)b;
  args.c = (வெறுமை *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (வெறுமை *)&alpha;
#இல்லையெனில்
  args.alpha = (வெறுமை *)alpha;
#endif

#எனில் !defined(COMPLEX) || defined(HEMM)
  args.beta  = (வெறுமை *)&beta;
#இல்லையெனில்
  args.beta  = (வெறுமை *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (Uplo == CblasUpper) uplo  = 0;
    எனில் (Uplo == CblasLower) uplo  = 1;

    எனில் (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    எனில் (Trans == CblasTrans)       trans = 1;
    எனில் (Trans == CblasConjNoTrans) trans = 0;
    எனில் (Trans == CblasConjTrans)   trans = 1;
#அல்லதெனில் !defined(HEMM)
    எனில் (Trans == CblasTrans)       trans = 1;
#இல்லையெனில்
    எனில் (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    எனில் (trans & 1) nrowa = args.k;

    எனில் (args.ldc < MAX(1,args.n)) info = 12;
    எனில் (args.ldb < MAX(1,nrowa))  info =  9;
    எனில் (args.lda < MAX(1,nrowa))  info =  7;
    எனில் (args.k < 0)               info =  4;
    எனில் (args.n < 0)               info =  3;
    எனில் (trans < 0)                info =  2;
    எனில் (uplo  < 0)                info =  1;
  }

  எனில் (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (வெறுமை *)CAlpha;
#endif

    எனில் (Uplo == CblasUpper) uplo  = 1;
    எனில் (Uplo == CblasLower) uplo  = 0;

    எனில் (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    எனில் (Trans == CblasTrans)       trans = 0;
    எனில் (Trans == CblasConjNoTrans) trans = 1;
    எனில் (Trans == CblasConjTrans)   trans = 0;
#அல்லதெனில் !defined(HEMM)
    எனில் (Trans == CblasTrans)       trans = 0;
#இல்லையெனில்
    எனில் (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    எனில் (trans & 1) nrowa = args.k;

    எனில் (args.ldc < MAX(1,args.n)) info = 12;
    எனில் (args.ldb < MAX(1,nrowa))  info =  9;
    எனில் (args.lda < MAX(1,nrowa))  info =  7;
    எனில் (args.k < 0)               info =  4;
    எனில் (args.n < 0)               info =  3;
    எனில் (trans < 0)                info =  2;
    எனில் (uplo  < 0)                info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

  எனில் (args.n == 0) திரும்பு;

#ifdef DYNAMIC_ARCH
வெளி எழுத்து* gotoblas_corename(வெறுமை);
#endif

#எனில் !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#எனில் defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
எனில் (strcmp(gotoblas_corename(), "armv9sme") == 0
#எனில் defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   எனில் (order == CblasRowMajor && n == ldc) {
     எனில் (Trans == CblasNoTrans && k == lda && k == ldb) {
      எனில் (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        திரும்பு;
      }இல்லையெனில் எனில் (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        திரும்பு;
      }
     }
     இல்லையெனில் எனில் (Trans == CblasTrans && n == lda && n ==ldb) {
      எனில் (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        திரும்பு;
      }இல்லையெனில் எனில் (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        திரும்பு;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  எனில் (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  எனில் (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } இல்லையெனில் {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = வெற்று;
  எனில் (args.n*args.k <1000)
	  args.nthreads =1 ;
  இல்லையெனில்
  args.nthreads = இலக்கம்_cpu_avail(3);

  எனில் (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, வெற்று, வெற்று, sa, sb, 0);

#ifdef SMP

  } இல்லையெனில் {

    syrk_thread(mode, &args, வெற்று, வெற்று, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  திரும்பு;
}
