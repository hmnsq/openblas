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
#မဟုတ်ပါက
#define ERROR_NAME "SSYR2K"
#endif
#မဟုတ်ပါက
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#မဟုတ်ပါက
#define ERROR_NAME "CSYR2K"
#endif
#မဟုတ်ပါက
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#မဟုတ်ပါက
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

static ကိန်းပြည့် (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#မဟုတ်ပါက
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ကိန်းပြည့် mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  ကိန်းပြည့် mode  =  BLAS_DOUBLE  | BLAS_REAL;
#မဟုတ်ပါက
  ကိန်းပြည့် mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#မဟုတ်ပါက
#ifdef XDOUBLE
  ကိန်းပြည့် mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  ကိန်းပြည့် mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#မဟုတ်ပါက
  ကိန်းပြည့် mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  ကိန်းပြည့် uplo;
  ကိန်းပြည့် trans;
  ကိန်းပြည့် nrowa;

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

  အကယ်၍ (uplo_arg  == 'U') uplo  = 0;
  အကယ်၍ (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  အကယ်၍ (trans_arg == 'N') trans = 0;
  အကယ်၍ (trans_arg == 'T') trans = 1;
  အကယ်၍ (trans_arg == 'C') trans = 1;
#မဟုတ်ပါက
#ifdef HEMM
  အကယ်၍ (trans_arg == 'N') trans = 0;
  အကယ်၍ (trans_arg == 'C') trans = 1;
#မဟုတ်ပါက
  အကယ်၍ (trans_arg == 'N') trans = 0;
  အကယ်၍ (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  အကယ်၍ (trans & 1) nrowa = args.k;

  info = 0;

  အကယ်၍ (args.ldc < MAX(1,args.n)) info = 12;
  အကယ်၍ (args.ldb < MAX(1,nrowa))  info =  9;
  အကယ်၍ (args.lda < MAX(1,nrowa))  info =  7;
  အကယ်၍ (args.k < 0)               info =  4;
  အကယ်၍ (args.n < 0)               info =  3;
  အကယ်၍ (trans < 0)                info =  2;
  အကယ်၍ (uplo  < 0)                info =  1;

  အကယ်၍ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#မဟုတ်ပါက

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, ga_nanaa CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#မဟုတ်ပါက
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
#endif
#အကယ်၍ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#မဟုတ်ပါက
	   void *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#မဟုတ်ပါက
	   void *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#အကယ်၍ !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ကိန်းပြည့် uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ကိန်းပြည့် mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  ကိန်းပြည့် mode  =  BLAS_DOUBLE  | BLAS_REAL;
#မဟုတ်ပါက
  ကိန်းပြည့် mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#မဟုတ်ပါက
#ifdef XDOUBLE
  ကိန်းပြည့် mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  ကိန်းပြည့် mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#မဟုတ်ပါက
  ကိန်းပြည့် mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
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
#မဟုတ်ပါက
  args.alpha = (void *)alpha;
#endif

#အကယ်၍ !defined(COMPLEX) || defined(HEMM)
  args.beta  = (void *)&beta;
#မဟုတ်ပါက
  args.beta  = (void *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  အကယ်၍ (order == CblasColMajor) {
    အကယ်၍ (Uplo == CblasUpper) uplo  = 0;
    အကယ်၍ (Uplo == CblasLower) uplo  = 1;

    အကယ်၍ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    အကယ်၍ (Trans == CblasTrans)       trans = 1;
    အကယ်၍ (Trans == CblasConjNoTrans) trans = 0;
    အကယ်၍ (Trans == CblasConjTrans)   trans = 1;
#elif !defined(HEMM)
    အကယ်၍ (Trans == CblasTrans)       trans = 1;
#မဟုတ်ပါက
    အကယ်၍ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    အကယ်၍ (trans & 1) nrowa = args.k;

    အကယ်၍ (args.ldc < MAX(1,args.n)) info = 12;
    အကယ်၍ (args.ldb < MAX(1,nrowa))  info =  9;
    အကယ်၍ (args.lda < MAX(1,nrowa))  info =  7;
    အကယ်၍ (args.k < 0)               info =  4;
    အကယ်၍ (args.n < 0)               info =  3;
    အကယ်၍ (trans < 0)                info =  2;
    အကယ်၍ (uplo  < 0)                info =  1;
  }

  အကယ်၍ (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (void *)CAlpha;
#endif

    အကယ်၍ (Uplo == CblasUpper) uplo  = 1;
    အကယ်၍ (Uplo == CblasLower) uplo  = 0;

    အကယ်၍ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    အကယ်၍ (Trans == CblasTrans)       trans = 0;
    အကယ်၍ (Trans == CblasConjNoTrans) trans = 1;
    အကယ်၍ (Trans == CblasConjTrans)   trans = 0;
#elif !defined(HEMM)
    အကယ်၍ (Trans == CblasTrans)       trans = 0;
#မဟုတ်ပါက
    အကယ်၍ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    အကယ်၍ (trans & 1) nrowa = args.k;

    အကယ်၍ (args.ldc < MAX(1,args.n)) info = 12;
    အကယ်၍ (args.ldb < MAX(1,nrowa))  info =  9;
    အကယ်၍ (args.lda < MAX(1,nrowa))  info =  7;
    အကယ်၍ (args.k < 0)               info =  4;
    အကယ်၍ (args.n < 0)               info =  3;
    အကယ်၍ (trans < 0)                info =  2;
    အကယ်၍ (uplo  < 0)                info =  1;
  }

  အကယ်၍ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

  အကယ်၍ (args.n == 0) ပြန်ပေး;

#ifdef DYNAMIC_ARCH
extern var_na* gotoblas_corename(void);
#endif

#အကယ်၍ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#အကယ်၍ defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#အကယ်၍ defined(DYNAMIC_ARCH)
အကယ်၍ (strcmp(gotoblas_corename(), "armv9sme") == 0
#အကယ်၍ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   အကယ်၍ (order == CblasRowMajor && n == ldc) {
     အကယ်၍ (Trans == CblasNoTrans && k == lda && k == ldb) {
      အကယ်၍ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ပြန်ပေး;
      }မဟုတ်ပါက အကယ်၍ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ပြန်ပေး;
      }
     }
     မဟုတ်ပါက အကယ်၍ (Trans == CblasTrans && n == lda && n ==ldb) {
      အကယ်၍ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ပြန်ပေး;
      }မဟုတ်ပါက အကယ်၍ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ပြန်ပေး;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  အကယ်၍ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  အကယ်၍ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } မဟုတ်ပါက {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = NULL;
  အကယ်၍ (args.n*args.k <1000)
	  args.nthreads =1 ;
  မဟုတ်ပါက
  args.nthreads = num_cpu_avail(3);

  အကယ်၍ (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } မဟုတ်ပါက {

    syrk_thread(mode, &args, NULL, NULL, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  ပြန်ပေး;
}
