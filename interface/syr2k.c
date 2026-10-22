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
#아니면
#define ERROR_NAME "SSYR2K"
#endif
#아니면
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#아니면
#define ERROR_NAME "CSYR2K"
#endif
#아니면
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#아니면
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

static 정수 (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#아니면
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
  정수 mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  정수 mode  =  BLAS_DOUBLE  | BLAS_REAL;
#아니면
  정수 mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#아니면
#ifdef XDOUBLE
  정수 mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  정수 mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#아니면
  정수 mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  정수 uplo;
  정수 trans;
  정수 nrowa;

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

  만약 (uplo_arg  == 'U') uplo  = 0;
  만약 (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  만약 (trans_arg == 'N') trans = 0;
  만약 (trans_arg == 'T') trans = 1;
  만약 (trans_arg == 'C') trans = 1;
#아니면
#ifdef HEMM
  만약 (trans_arg == 'N') trans = 0;
  만약 (trans_arg == 'C') trans = 1;
#아니면
  만약 (trans_arg == 'N') trans = 0;
  만약 (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  만약 (trans & 1) nrowa = args.k;

  info = 0;

  만약 (args.ldc < MAX(1,args.n)) info = 12;
  만약 (args.ldb < MAX(1,nrowa))  info =  9;
  만약 (args.lda < MAX(1,nrowa))  info =  7;
  만약 (args.k < 0)               info =  4;
  만약 (args.n < 0)               info =  3;
  만약 (trans < 0)                info =  2;
  만약 (uplo  < 0)                info =  1;

  만약 (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    반환;
  }

#아니면

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#아니면
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
#endif
#만약 !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#아니면
	   void *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#아니면
	   void *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#만약 !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  정수 uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  정수 mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  정수 mode  =  BLAS_DOUBLE  | BLAS_REAL;
#아니면
  정수 mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#아니면
#ifdef XDOUBLE
  정수 mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  정수 mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#아니면
  정수 mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
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
#아니면
  args.alpha = (void *)alpha;
#endif

#만약 !defined(COMPLEX) || defined(HEMM)
  args.beta  = (void *)&beta;
#아니면
  args.beta  = (void *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  만약 (order == CblasColMajor) {
    만약 (Uplo == CblasUpper) uplo  = 0;
    만약 (Uplo == CblasLower) uplo  = 1;

    만약 (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    만약 (Trans == CblasTrans)       trans = 1;
    만약 (Trans == CblasConjNoTrans) trans = 0;
    만약 (Trans == CblasConjTrans)   trans = 1;
#elif !defined(HEMM)
    만약 (Trans == CblasTrans)       trans = 1;
#아니면
    만약 (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    만약 (trans & 1) nrowa = args.k;

    만약 (args.ldc < MAX(1,args.n)) info = 12;
    만약 (args.ldb < MAX(1,nrowa))  info =  9;
    만약 (args.lda < MAX(1,nrowa))  info =  7;
    만약 (args.k < 0)               info =  4;
    만약 (args.n < 0)               info =  3;
    만약 (trans < 0)                info =  2;
    만약 (uplo  < 0)                info =  1;
  }

  만약 (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (void *)CAlpha;
#endif

    만약 (Uplo == CblasUpper) uplo  = 1;
    만약 (Uplo == CblasLower) uplo  = 0;

    만약 (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    만약 (Trans == CblasTrans)       trans = 0;
    만약 (Trans == CblasConjNoTrans) trans = 1;
    만약 (Trans == CblasConjTrans)   trans = 0;
#elif !defined(HEMM)
    만약 (Trans == CblasTrans)       trans = 0;
#아니면
    만약 (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    만약 (trans & 1) nrowa = args.k;

    만약 (args.ldc < MAX(1,args.n)) info = 12;
    만약 (args.ldb < MAX(1,nrowa))  info =  9;
    만약 (args.lda < MAX(1,nrowa))  info =  7;
    만약 (args.k < 0)               info =  4;
    만약 (args.n < 0)               info =  3;
    만약 (trans < 0)                info =  2;
    만약 (uplo  < 0)                info =  1;
  }

  만약 (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    반환;
  }

  만약 (args.n == 0) 반환;

#ifdef DYNAMIC_ARCH
extern char* gotoblas_corename(void);
#endif

#만약 !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#만약 defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#만약 defined(DYNAMIC_ARCH)
만약 (strcmp(gotoblas_corename(), "armv9sme") == 0
#만약 defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   만약 (order == CblasRowMajor && n == ldc) {
     만약 (Trans == CblasNoTrans && k == lda && k == ldb) {
      만약 (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        반환;
      }아니면 만약 (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        반환;
      }
     }
     아니면 만약 (Trans == CblasTrans && n == lda && n ==ldb) {
      만약 (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        반환;
      }아니면 만약 (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        반환;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  만약 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    반환;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  만약 (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } 아니면 {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = NULL;
  만약 (args.n*args.k <1000)
	  args.nthreads =1 ;
  아니면
  args.nthreads = num_cpu_avail(3);

  만약 (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } 아니면 {

    syrk_thread(mode, &args, NULL, NULL, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  반환;
}
