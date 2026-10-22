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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#elif defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#инако
#define ERROR_NAME "SSYR2K"
#endif
#инако
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#инако
#define ERROR_NAME "CSYR2K"
#endif
#инако
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#инако
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

static цео (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#инако
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, char *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  char uplo_arg  = *UPLO;
  char trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  цео mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  цео mode  =  BLAS_DOUBLE  | BLAS_REAL;
#инако
  цео mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#инако
#ifdef XDOUBLE
  цео mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  цео mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#инако
  цео mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  цео uplo;
  цео trans;
  цео nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (void *)a;
  args.b = (void *)b;
  args.c = (void *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (void *)alpha;
  args.beta  = (void *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ако (uplo_arg  == 'U') uplo  = 0;
  ако (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  ако (trans_arg == 'N') trans = 0;
  ако (trans_arg == 'T') trans = 1;
  ако (trans_arg == 'C') trans = 1;
#инако
#ifdef HEMM
  ако (trans_arg == 'N') trans = 0;
  ако (trans_arg == 'C') trans = 1;
#инако
  ако (trans_arg == 'N') trans = 0;
  ако (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  ако (trans & 1) nrowa = args.k;

  info = 0;

  ако (args.ldc < MAX(1,args.n)) info = 12;
  ако (args.ldb < MAX(1,nrowa))  info =  9;
  ако (args.lda < MAX(1,nrowa))  info =  7;
  ако (args.k < 0)               info =  4;
  ако (args.n < 0)               info =  3;
  ако (trans < 0)                info =  2;
  ако (uplo  < 0)                info =  1;

  ако (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    врати;
  }

#инако

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#инако
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
#endif
#ако !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#инако
	   void *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#инако
	   void *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#ако !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  цео uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  цео mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  цео mode  =  BLAS_DOUBLE  | BLAS_REAL;
#инако
  цео mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#инако
#ifdef XDOUBLE
  цео mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  цео mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#инако
  цео mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (void *)a;
  args.b = (void *)b;
  args.c = (void *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (void *)&alpha;
#инако
  args.alpha = (void *)alpha;
#endif

#ако !defined(COMPLEX) || defined(HEMM)
  args.beta  = (void *)&beta;
#инако
  args.beta  = (void *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  ако (order == CblasColMajor) {
    ако (Uplo == CblasUpper) uplo  = 0;
    ако (Uplo == CblasLower) uplo  = 1;

    ако (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    ако (Trans == CblasTrans)       trans = 1;
    ако (Trans == CblasConjNoTrans) trans = 0;
    ако (Trans == CblasConjTrans)   trans = 1;
#elif !defined(HEMM)
    ако (Trans == CblasTrans)       trans = 1;
#инако
    ако (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    ако (trans & 1) nrowa = args.k;

    ако (args.ldc < MAX(1,args.n)) info = 12;
    ако (args.ldb < MAX(1,nrowa))  info =  9;
    ако (args.lda < MAX(1,nrowa))  info =  7;
    ако (args.k < 0)               info =  4;
    ако (args.n < 0)               info =  3;
    ако (trans < 0)                info =  2;
    ако (uplo  < 0)                info =  1;
  }

  ако (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (void *)CAlpha;
#endif

    ако (Uplo == CblasUpper) uplo  = 1;
    ако (Uplo == CblasLower) uplo  = 0;

    ако (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    ако (Trans == CblasTrans)       trans = 0;
    ако (Trans == CblasConjNoTrans) trans = 1;
    ако (Trans == CblasConjTrans)   trans = 0;
#elif !defined(HEMM)
    ако (Trans == CblasTrans)       trans = 0;
#инако
    ако (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    ако (trans & 1) nrowa = args.k;

    ако (args.ldc < MAX(1,args.n)) info = 12;
    ако (args.ldb < MAX(1,nrowa))  info =  9;
    ако (args.lda < MAX(1,nrowa))  info =  7;
    ако (args.k < 0)               info =  4;
    ако (args.n < 0)               info =  3;
    ако (trans < 0)                info =  2;
    ако (uplo  < 0)                info =  1;
  }

  ако (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    врати;
  }

  ако (args.n == 0) врати;

#ifdef DYNAMIC_ARCH
extern char* gotoblas_corename(void);
#endif

#ако !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ако defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ако defined(DYNAMIC_ARCH)
ако (strcmp(gotoblas_corename(), "armv9sme") == 0
#ако defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ако (order == CblasRowMajor && n == ldc) {
     ако (Trans == CblasNoTrans && k == lda && k == ldb) {
      ако (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        врати;
      }инако ако (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        врати;
      }
     }
     инако ако (Trans == CblasTrans && n == lda && n ==ldb) {
      ако (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        врати;
      }инако ако (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        врати;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ако (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    врати;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  ако (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } инако {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = NULL;
  ако (args.n*args.k <1000)
	  args.nthreads =1 ;
  инако
  args.nthreads = num_cpu_avail(3);

  ако (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } инако {

    syrk_thread(mode, &args, NULL, NULL, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  врати;
}
