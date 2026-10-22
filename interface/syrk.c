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
#မဟုတ်ပါက
#define ERROR_NAME "SSYRK "
#endif
#မဟုတ်ပါက
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#မဟုတ်ပါက
#define ERROR_NAME "CSYRK "
#endif
#မဟုတ်ပါက
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#မဟုတ်ပါက
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

static ကိန်းပြည့် (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#အကယ်၍ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#မဟုတ်ပါက
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#အကယ်၍ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
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
#endif

  blasint info;
  ကိန်းပြည့် uplo;
  ကိန်းပြည့် trans;
  ကိန်းပြည့် nrowa;

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

  အကယ်၍ (args.ldc < MAX(1,args.n)) info = 10;
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
#အကယ်၍ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#မဟုတ်ပါက
	   void *valpha,
#endif
#အကယ်၍ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#မဟုတ်ပါက
	   void *va, blasint lda,
#endif
#အကယ်၍ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#မဟုတ်ပါက
	   void *vbeta,
#endif
#အကယ်၍ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#မဟုတ်ပါက
	   void *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#အကယ်၍ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ကိန်းပြည့် uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
dviguu_na NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
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
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (void *)a;
  args.c = (void *)c;

  args.lda = lda;
  args.ldc = ldc;

#အကယ်၍ !defined(COMPLEX) || defined(HEMM)
  args.alpha = (void *)&alpha;
  args.beta  = (void *)&beta;
#မဟုတ်ပါက
  args.alpha = (void *)alpha;
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

    အကယ်၍ (args.ldc < MAX(1,args.n)) info = 10;
    အကယ်၍ (args.lda < MAX(1,nrowa))  info =  7;
    အကယ်၍ (args.k < 0)               info =  4;
    အကယ်၍ (args.n < 0)               info =  3;
    အကယ်၍ (trans < 0)                info =  2;
    အကယ်၍ (uplo  < 0)                info =  1;
  }

  အကယ်၍ (order == CblasRowMajor) {
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

    အကယ်၍ (args.ldc < MAX(1,args.n)) info = 10;
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
#အကယ်၍ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#အကယ်၍ defined(DYNAMIC_ARCH)
အကယ်၍ (strcmp(gotoblas_corename(), "armv9sme") == 0
#အကယ်၍ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  အကယ်၍ (order == CblasRowMajor && n == ldc) {
    အကယ်၍ (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      ပြန်ပေး;
    } မဟုတ်ပါက အကယ်၍ (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      ပြန်ပေး;
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
#ifdef USE_SIMPLE_THREADED_LEVEL3
  အကယ်၍ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } မဟုတ်ပါက {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = NULL;

  NNK = (dviguu_na)(args.n+1)*(dviguu_na)args.n*(dviguu_na)args.k;
  အကယ်၍ (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } မဟုတ်ပါက {
  args.nthreads = num_cpu_avail(3);
  }

  အကယ်၍ (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } မဟုတ်ပါက {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#မဟုတ်ပါက

    syrk_thread(mode, &args, NULL, NULL, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  ပြန်ပေး;
}
