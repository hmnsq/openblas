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
#དེ་མིན
#define ERROR_NAME "SSYR2K"
#endif
#དེ་མིན
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#དེ་མིན
#define ERROR_NAME "CSYR2K"
#endif
#དེ་མིན
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#དེ་མིན
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

static ཧྲིལ (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#དེ་མིན
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
  ཧྲིལ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  ཧྲིལ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#དེ་མིན
  ཧྲིལ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#དེ་མིན
#ifdef XDOUBLE
  ཧྲིལ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  ཧྲིལ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#དེ་མིན
  ཧྲིལ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  ཧྲིལ uplo;
  ཧྲིལ trans;
  ཧྲིལ nrowa;

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

  གལ་ཏེ (uplo_arg  == 'U') uplo  = 0;
  གལ་ཏེ (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  གལ་ཏེ (trans_arg == 'N') trans = 0;
  གལ་ཏེ (trans_arg == 'T') trans = 1;
  གལ་ཏེ (trans_arg == 'C') trans = 1;
#དེ་མིན
#ifdef HEMM
  གལ་ཏེ (trans_arg == 'N') trans = 0;
  གལ་ཏེ (trans_arg == 'C') trans = 1;
#དེ་མིན
  གལ་ཏེ (trans_arg == 'N') trans = 0;
  གལ་ཏེ (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  གལ་ཏེ (trans & 1) nrowa = args.k;

  info = 0;

  གལ་ཏེ (args.ldc < MAX(1,args.n)) info = 12;
  གལ་ཏེ (args.ldb < MAX(1,nrowa))  info =  9;
  གལ་ཏེ (args.lda < MAX(1,nrowa))  info =  7;
  གལ་ཏེ (args.k < 0)               info =  4;
  གལ་ཏེ (args.n < 0)               info =  3;
  གལ་ཏེ (trans < 0)                info =  2;
  གལ་ཏེ (uplo  < 0)                info =  1;

  གལ་ཏེ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ཕྱིར;
  }

#དེ་མིན

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, ga_nanaa CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#དེ་མིན
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
#endif
#གལ་ཏེ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#དེ་མིན
	   void *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#དེ་མིན
	   void *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#གལ་ཏེ !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ཧྲིལ uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ཧྲིལ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  ཧྲིལ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#དེ་མིན
  ཧྲིལ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#དེ་མིན
#ifdef XDOUBLE
  ཧྲིལ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  ཧྲིལ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#དེ་མིན
  ཧྲིལ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
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
#དེ་མིན
  args.alpha = (void *)alpha;
#endif

#གལ་ཏེ !defined(COMPLEX) || defined(HEMM)
  args.beta  = (void *)&beta;
#དེ་མིན
  args.beta  = (void *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  གལ་ཏེ (order == CblasColMajor) {
    གལ་ཏེ (Uplo == CblasUpper) uplo  = 0;
    གལ་ཏེ (Uplo == CblasLower) uplo  = 1;

    གལ་ཏེ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    གལ་ཏེ (Trans == CblasTrans)       trans = 1;
    གལ་ཏེ (Trans == CblasConjNoTrans) trans = 0;
    གལ་ཏེ (Trans == CblasConjTrans)   trans = 1;
#elif !defined(HEMM)
    གལ་ཏེ (Trans == CblasTrans)       trans = 1;
#དེ་མིན
    གལ་ཏེ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    གལ་ཏེ (trans & 1) nrowa = args.k;

    གལ་ཏེ (args.ldc < MAX(1,args.n)) info = 12;
    གལ་ཏེ (args.ldb < MAX(1,nrowa))  info =  9;
    གལ་ཏེ (args.lda < MAX(1,nrowa))  info =  7;
    གལ་ཏེ (args.k < 0)               info =  4;
    གལ་ཏེ (args.n < 0)               info =  3;
    གལ་ཏེ (trans < 0)                info =  2;
    གལ་ཏེ (uplo  < 0)                info =  1;
  }

  གལ་ཏེ (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (void *)CAlpha;
#endif

    གལ་ཏེ (Uplo == CblasUpper) uplo  = 1;
    གལ་ཏེ (Uplo == CblasLower) uplo  = 0;

    གལ་ཏེ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    གལ་ཏེ (Trans == CblasTrans)       trans = 0;
    གལ་ཏེ (Trans == CblasConjNoTrans) trans = 1;
    གལ་ཏེ (Trans == CblasConjTrans)   trans = 0;
#elif !defined(HEMM)
    གལ་ཏེ (Trans == CblasTrans)       trans = 0;
#དེ་མིན
    གལ་ཏེ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    གལ་ཏེ (trans & 1) nrowa = args.k;

    གལ་ཏེ (args.ldc < MAX(1,args.n)) info = 12;
    གལ་ཏེ (args.ldb < MAX(1,nrowa))  info =  9;
    གལ་ཏེ (args.lda < MAX(1,nrowa))  info =  7;
    གལ་ཏེ (args.k < 0)               info =  4;
    གལ་ཏེ (args.n < 0)               info =  3;
    གལ་ཏེ (trans < 0)                info =  2;
    གལ་ཏེ (uplo  < 0)                info =  1;
  }

  གལ་ཏེ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ཕྱིར;
  }

  གལ་ཏེ (args.n == 0) ཕྱིར;

#ifdef DYNAMIC_ARCH
extern var_na* gotoblas_corename(void);
#endif

#གལ་ཏེ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#གལ་ཏེ defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#གལ་ཏེ defined(DYNAMIC_ARCH)
གལ་ཏེ (strcmp(gotoblas_corename(), "armv9sme") == 0
#གལ་ཏེ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   གལ་ཏེ (order == CblasRowMajor && n == ldc) {
     གལ་ཏེ (Trans == CblasNoTrans && k == lda && k == ldb) {
      གལ་ཏེ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ཕྱིར;
      }དེ་མིན གལ་ཏེ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ཕྱིར;
      }
     }
     དེ་མིན གལ་ཏེ (Trans == CblasTrans && n == lda && n ==ldb) {
      གལ་ཏེ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ཕྱིར;
      }དེ་མིན གལ་ཏེ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ཕྱིར;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  གལ་ཏེ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ཕྱིར;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  གལ་ཏེ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } དེ་མིན {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = NULL;
  གལ་ཏེ (args.n*args.k <1000)
	  args.nthreads =1 ;
  དེ་མིན
  args.nthreads = num_cpu_avail(3);

  གལ་ཏེ (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } དེ་མིན {

    syrk_thread(mode, &args, NULL, NULL, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  ཕྱིར;
}
