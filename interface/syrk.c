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
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#elif defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#དེ་མིན
#define ERROR_NAME "SSYRK "
#endif
#དེ་མིན
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#དེ་མིན
#define ERROR_NAME "CSYRK "
#endif
#དེ་མིན
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#དེ་མིན
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

static ཧྲིལ (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#གལ་ཏེ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#དེ་མིན
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#གལ་ཏེ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  dviguu_na NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
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
#endif

  blasint info;
  ཧྲིལ uplo;
  ཧྲིལ trans;
  ཧྲིལ nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (void *)a;
  args.c = (void *)c;

  args.lda = *ldA;
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

  གལ་ཏེ (args.ldc < MAX(1,args.n)) info = 10;
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
#གལ་ཏེ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#དེ་མིན
	   void *valpha,
#endif
#གལ་ཏེ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#དེ་མིན
	   void *va, blasint lda,
#endif
#གལ་ཏེ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#དེ་མིན
	   void *vbeta,
#endif
#གལ་ཏེ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#དེ་མིན
	   void *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#གལ་ཏེ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ཧྲིལ uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
dviguu_na NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
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
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (void *)a;
  args.c = (void *)c;

  args.lda = lda;
  args.ldc = ldc;

#གལ་ཏེ !defined(COMPLEX) || defined(HEMM)
  args.alpha = (void *)&alpha;
  args.beta  = (void *)&beta;
#དེ་མིན
  args.alpha = (void *)alpha;
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

    གལ་ཏེ (args.ldc < MAX(1,args.n)) info = 10;
    གལ་ཏེ (args.lda < MAX(1,nrowa))  info =  7;
    གལ་ཏེ (args.k < 0)               info =  4;
    གལ་ཏེ (args.n < 0)               info =  3;
    གལ་ཏེ (trans < 0)                info =  2;
    གལ་ཏེ (uplo  < 0)                info =  1;
  }

  གལ་ཏེ (order == CblasRowMajor) {
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

    གལ་ཏེ (args.ldc < MAX(1,args.n)) info = 10;
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
#གལ་ཏེ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#གལ་ཏེ defined(DYNAMIC_ARCH)
གལ་ཏེ (strcmp(gotoblas_corename(), "armv9sme") == 0
#གལ་ཏེ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  གལ་ཏེ (order == CblasRowMajor && n == ldc) {
    གལ་ཏེ (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      ཕྱིར;
    } དེ་མིན གལ་ཏེ (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      ཕྱིར;
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
#ifdef USE_SIMPLE_THREADED_LEVEL3
  གལ་ཏེ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } དེ་མིན {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = NULL;

  NNK = (dviguu_na)(args.n+1)*(dviguu_na)args.n*(dviguu_na)args.k;
  གལ་ཏེ (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } དེ་མིན {
  args.nthreads = num_cpu_avail(3);
  }

  གལ་ཏེ (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } དེ་མིན {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#དེ་མིན

    syrk_thread(mode, &args, NULL, NULL, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  ཕྱིར;
}
