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
#иначе
#define ERROR_NAME "SSYRK "
#endif
#иначе
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#иначе
#define ERROR_NAME "CSYRK "
#endif
#иначе
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#иначе
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

static целое (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#если defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#иначе
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#если defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, char *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  char uplo_arg  = *UPLO;
  char trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  double NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  целое mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  целое mode  =  BLAS_DOUBLE  | BLAS_REAL;
#иначе
  целое mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#иначе
#ifdef XDOUBLE
  целое mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  целое mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#иначе
  целое mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  целое uplo;
  целое trans;
  целое nrowa;

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

  если (uplo_arg  == 'U') uplo  = 0;
  если (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  если (trans_arg == 'N') trans = 0;
  если (trans_arg == 'T') trans = 1;
  если (trans_arg == 'C') trans = 1;
#иначе
#ifdef HEMM
  если (trans_arg == 'N') trans = 0;
  если (trans_arg == 'C') trans = 1;
#иначе
  если (trans_arg == 'N') trans = 0;
  если (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  если (trans & 1) nrowa = args.k;

  info = 0;

  если (args.ldc < MAX(1,args.n)) info = 10;
  если (args.lda < MAX(1,nrowa))  info =  7;
  если (args.k < 0)               info =  4;
  если (args.n < 0)               info =  3;
  если (trans < 0)                info =  2;
  если (uplo  < 0)                info =  1;

  если (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#иначе

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#если !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#иначе
	   void *valpha,
#endif
#если !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#иначе
	   void *va, blasint lda,
#endif
#если !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#иначе
	   void *vbeta,
#endif
#если !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#иначе
	   void *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#если !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  целое uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
double NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  целое mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  целое mode  =  BLAS_DOUBLE  | BLAS_REAL;
#иначе
  целое mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#иначе
#ifdef XDOUBLE
  целое mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  целое mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#иначе
  целое mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
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

#если !defined(COMPLEX) || defined(HEMM)
  args.alpha = (void *)&alpha;
  args.beta  = (void *)&beta;
#иначе
  args.alpha = (void *)alpha;
  args.beta  = (void *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  если (order == CblasColMajor) {
    если (Uplo == CblasUpper) uplo  = 0;
    если (Uplo == CblasLower) uplo  = 1;

    если (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    если (Trans == CblasTrans)       trans = 1;
    если (Trans == CblasConjNoTrans) trans = 0;
    если (Trans == CblasConjTrans)   trans = 1;
#elif !defined(HEMM)
    если (Trans == CblasTrans)       trans = 1;
#иначе
    если (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    если (trans & 1) nrowa = args.k;

    если (args.ldc < MAX(1,args.n)) info = 10;
    если (args.lda < MAX(1,nrowa))  info =  7;
    если (args.k < 0)               info =  4;
    если (args.n < 0)               info =  3;
    если (trans < 0)                info =  2;
    если (uplo  < 0)                info =  1;
  }

  если (order == CblasRowMajor) {
    если (Uplo == CblasUpper) uplo  = 1;
    если (Uplo == CblasLower) uplo  = 0;

    если (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    если (Trans == CblasTrans)       trans = 0;
    если (Trans == CblasConjNoTrans) trans = 1;
    если (Trans == CblasConjTrans)   trans = 0;
#elif !defined(HEMM)
    если (Trans == CblasTrans)       trans = 0;
#иначе
    если (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    если (trans & 1) nrowa = args.k;

    если (args.ldc < MAX(1,args.n)) info = 10;
    если (args.lda < MAX(1,nrowa))  info =  7;
    если (args.k < 0)               info =  4;
    если (args.n < 0)               info =  3;
    если (trans < 0)                info =  2;
    если (uplo  < 0)                info =  1;
  }

  если (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

  если (args.n == 0) вернуть;

#ifdef DYNAMIC_ARCH
extern char* gotoblas_corename(void);
#endif
	
#если !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#если defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#если defined(DYNAMIC_ARCH)
если (strcmp(gotoblas_corename(), "armv9sme") == 0
#если defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  если (order == CblasRowMajor && n == ldc) {
    если (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      вернуть;
    } иначе если (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      вернуть;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  если (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  если (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } иначе {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = NULL;

  NNK = (double)(args.n+1)*(double)args.n*(double)args.k;
  если (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } иначе {
  args.nthreads = num_cpu_avail(3);
  }

  если (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } иначе {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, NULL, NULL, sa, sb, 0);

#иначе

    syrk_thread(mode, &args, NULL, NULL, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  вернуть;
}
