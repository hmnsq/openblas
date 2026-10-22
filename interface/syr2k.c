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

#شامل <stdio.h>
#شامل <ctype.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#وگرنه
#define ERROR_NAME "SSYR2K"
#endif
#وگرنه
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#وگرنه
#define ERROR_NAME "CSYR2K"
#endif
#وگرنه
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#وگرنه
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

ایستا صحیح (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#وگرنه
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *UPLO, نویسه *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  نویسه uplo_arg  = *UPLO;
  نویسه trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  صحیح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنهاگر defined(DOUBLE)
  صحیح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وگرنه
  صحیح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وگرنه
#ifdef XDOUBLE
  صحیح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنهاگر defined(DOUBLE)
  صحیح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وگرنه
  صحیح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  صحیح uplo;
  صحیح trans;
  صحیح nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (تهی *)a;
  args.b = (تهی *)b;
  args.c = (تهی *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (تهی *)alpha;
  args.beta  = (تهی *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
  اگر (trans_arg == 'C') trans = 1;
#وگرنه
#ifdef HEMM
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'C') trans = 1;
#وگرنه
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  اگر (trans & 1) nrowa = args.k;

  info = 0;

  اگر (args.ldc < MAX(1,args.n)) info = 12;
  اگر (args.ldb < MAX(1,nrowa))  info =  9;
  اگر (args.lda < MAX(1,nrowa))  info =  7;
  اگر (args.k < 0)               info =  4;
  اگر (args.n < 0)               info =  3;
  اگر (trans < 0)                info =  2;
  اگر (uplo  < 0)                info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order, شمارش CBLAS_UPLO Uplo, شمارش CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#وگرنه
	   تهی *valpha,
	   تهی *va, blasint lda,
	   تهی *vb, blasint ldb,
#endif
#اگر !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وگرنه
	   تهی *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#وگرنه
	   تهی *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#اگر !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  صحیح uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  صحیح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنهاگر defined(DOUBLE)
  صحیح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وگرنه
  صحیح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وگرنه
#ifdef XDOUBLE
  صحیح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنهاگر defined(DOUBLE)
  صحیح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وگرنه
  صحیح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (تهی *)a;
  args.b = (تهی *)b;
  args.c = (تهی *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (تهی *)&alpha;
#وگرنه
  args.alpha = (تهی *)alpha;
#endif

#اگر !defined(COMPLEX) || defined(HEMM)
  args.beta  = (تهی *)&beta;
#وگرنه
  args.beta  = (تهی *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  اگر (order == CblasColMajor) {
    اگر (Uplo == CblasUpper) uplo  = 0;
    اگر (Uplo == CblasLower) uplo  = 1;

    اگر (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    اگر (Trans == CblasTrans)       trans = 1;
    اگر (Trans == CblasConjNoTrans) trans = 0;
    اگر (Trans == CblasConjTrans)   trans = 1;
#وگرنهاگر !defined(HEMM)
    اگر (Trans == CblasTrans)       trans = 1;
#وگرنه
    اگر (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    اگر (trans & 1) nrowa = args.k;

    اگر (args.ldc < MAX(1,args.n)) info = 12;
    اگر (args.ldb < MAX(1,nrowa))  info =  9;
    اگر (args.lda < MAX(1,nrowa))  info =  7;
    اگر (args.k < 0)               info =  4;
    اگر (args.n < 0)               info =  3;
    اگر (trans < 0)                info =  2;
    اگر (uplo  < 0)                info =  1;
  }

  اگر (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (تهی *)CAlpha;
#endif

    اگر (Uplo == CblasUpper) uplo  = 1;
    اگر (Uplo == CblasLower) uplo  = 0;

    اگر (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    اگر (Trans == CblasTrans)       trans = 0;
    اگر (Trans == CblasConjNoTrans) trans = 1;
    اگر (Trans == CblasConjTrans)   trans = 0;
#وگرنهاگر !defined(HEMM)
    اگر (Trans == CblasTrans)       trans = 0;
#وگرنه
    اگر (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    اگر (trans & 1) nrowa = args.k;

    اگر (args.ldc < MAX(1,args.n)) info = 12;
    اگر (args.ldb < MAX(1,nrowa))  info =  9;
    اگر (args.lda < MAX(1,nrowa))  info =  7;
    اگر (args.k < 0)               info =  4;
    اگر (args.n < 0)               info =  3;
    اگر (trans < 0)                info =  2;
    اگر (uplo  < 0)                info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

  اگر (args.n == 0) بازگشت;

#ifdef DYNAMIC_ARCH
بیرونی نویسه* gotoblas_corename(تهی);
#endif

#اگر !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#اگر defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
اگر (strcmp(gotoblas_corename(), "armv9sme") == 0
#اگر defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   اگر (order == CblasRowMajor && n == ldc) {
     اگر (Trans == CblasNoTrans && k == lda && k == ldb) {
      اگر (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        بازگشت;
      }وگرنه اگر (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        بازگشت;
      }
     }
     وگرنه اگر (Trans == CblasTrans && n == lda && n ==ldb) {
      اگر (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        بازگشت;
      }وگرنه اگر (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        بازگشت;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  اگر (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  اگر (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وگرنه {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = پوچ;
  اگر (args.n*args.k <1000)
	  args.nthreads =1 ;
  وگرنه
  args.nthreads = رقم_cpu_avail(3);

  اگر (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, پوچ, پوچ, sa, sb, 0);

#ifdef SMP

  } وگرنه {

    syrk_thread(mode, &args, پوچ, پوچ, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  بازگشت;
}
