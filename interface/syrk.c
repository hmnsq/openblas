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
#else
#define ERROR_NAME "SSYRK "
#endif
#else
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#else
#define ERROR_NAME "CSYRK "
#endif
#else
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#else
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

static int (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#else
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
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

  blas_arg_t मगा_कातविरान;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  double NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  int mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  int mode  =  BLAS_DOUBLE  | BLAS_REAL;
#else
  int mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#else
#ifdef XDOUBLE
  int mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  int mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#else
  int mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  int uplo;
  int trans;
  int nrowa;

  इलिम्बाग_DEBUG_NAME;

  मगा_कातविरान.n = *N;
  मगा_कातविरान.k = *K;

  मगा_कातविरान.a = (void *)a;
  मगा_कातविरान.c = (void *)c;

  मगा_कातविरान.lda = *ldA;
  मगा_कातविरान.ldc = *ldC;

  मगा_कातविरान.alpha = (void *)alpha;
  मगा_कातविरान.beta  = (void *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  if (uplo_arg  == 'U') uplo  = 0;
  if (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  if (trans_arg == 'N') trans = 0;
  if (trans_arg == 'T') trans = 1;
  if (trans_arg == 'C') trans = 1;
#else
#ifdef HEMM
  if (trans_arg == 'N') trans = 0;
  if (trans_arg == 'C') trans = 1;
#else
  if (trans_arg == 'N') trans = 0;
  if (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = मगा_कातविरान.n;
  if (trans & 1) nrowa = मगा_कातविरान.k;

  info = 0;

  if (मगा_कातविरान.ldc < MAX(1,मगा_कातविरान.n)) info = 10;
  if (मगा_कातविरान.lda < MAX(1,nrowa))  info =  7;
  if (मगा_कातविरान.k < 0)               info =  4;
  if (मगा_कातविरान.n < 0)               info =  3;
  if (trans < 0)                info =  2;
  if (uplo  < 0)                info =  1;

  if (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

#else

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#if !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#else
	   void *valpha,
#endif
#if !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#else
	   void *va, blasint lda,
#endif
#if !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#else
	   void *vbeta,
#endif
#if !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#else
	   void *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#if !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t मगा_कातविरान;
  int uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
double NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  int mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  int mode  =  BLAS_DOUBLE  | BLAS_REAL;
#else
  int mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#else
#ifdef XDOUBLE
  int mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  int mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#else
  int mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  इलिम्बाग_DEBUG_CNAME;

  मगा_कातविरान.n = n;
  मगा_कातविरान.k = k;

  मगा_कातविरान.a = (void *)a;
  मगा_कातविरान.c = (void *)c;

  मगा_कातविरान.lda = lda;
  मगा_कातविरान.ldc = ldc;

#if !defined(COMPLEX) || defined(HEMM)
  मगा_कातविरान.alpha = (void *)&alpha;
  मगा_कातविरान.beta  = (void *)&beta;
#else
  मगा_कातविरान.alpha = (void *)alpha;
  मगा_कातविरान.beta  = (void *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  if (order == CblasColMajor) {
    if (Uplo == CblasUpper) uplo  = 0;
    if (Uplo == CblasLower) uplo  = 1;

    if (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    if (Trans == CblasTrans)       trans = 1;
    if (Trans == CblasConjNoTrans) trans = 0;
    if (Trans == CblasConjTrans)   trans = 1;
#elif !defined(HEMM)
    if (Trans == CblasTrans)       trans = 1;
#else
    if (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = मगा_कातविरान.n;
    if (trans & 1) nrowa = मगा_कातविरान.k;

    if (मगा_कातविरान.ldc < MAX(1,मगा_कातविरान.n)) info = 10;
    if (मगा_कातविरान.lda < MAX(1,nrowa))  info =  7;
    if (मगा_कातविरान.k < 0)               info =  4;
    if (मगा_कातविरान.n < 0)               info =  3;
    if (trans < 0)                info =  2;
    if (uplo  < 0)                info =  1;
  }

  if (order == CblasRowMajor) {
    if (Uplo == CblasUpper) uplo  = 1;
    if (Uplo == CblasLower) uplo  = 0;

    if (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    if (Trans == CblasTrans)       trans = 0;
    if (Trans == CblasConjNoTrans) trans = 1;
    if (Trans == CblasConjTrans)   trans = 0;
#elif !defined(HEMM)
    if (Trans == CblasTrans)       trans = 0;
#else
    if (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = मगा_कातविरान.n;
    if (trans & 1) nrowa = मगा_कातविरान.k;

    if (मगा_कातविरान.ldc < MAX(1,मगा_कातविरान.n)) info = 10;
    if (मगा_कातविरान.lda < MAX(1,nrowa))  info =  7;
    if (मगा_कातविरान.k < 0)               info =  4;
    if (मगा_कातविरान.n < 0)               info =  3;
    if (trans < 0)                info =  2;
    if (uplo  < 0)                info =  1;
  }

  if (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

  if (मगा_कातविरान.n == 0) return;

#ifdef DYNAMIC_ARCH
extern char* gotoblas_corename(void);
#endif
	
#if !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#if defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#if defined(DYNAMIC_ARCH)
if (strcmp(gotoblas_corename(), "armv9sme") == 0
#if defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  if (order == CblasRowMajor && n == ldc) {
    if (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      return;
    } else if (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      return;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  if (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  if (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } else {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  मगा_कातविरान.common = NULL;

  NNK = (double)(मगा_कातविरान.n+1)*(double)मगा_कातविरान.n*(double)मगा_कातविरान.k;
  if (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  मगा_कातविरान.nthreads = 1;
  } else {
  मगा_कातविरान.nthreads = नुमेरो_cpu_avail(3);
  }

  if (मगा_कातविरान.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } else {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

#else

    syrk_thread(mode, &मगा_कातविरान, NULL, NULL, syrk[(uplo << 1) | trans ], sa, sb, मगा_कातविरान.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, मगा_कातविरान.n * मगा_कातविरान.k + मगा_कातविरान.n * मगा_कातविरान.n / 2, मगा_कातविरान.n * मगा_कातविरान.n * मगा_कातविरान.k);

  IDEBUG_END;

  return;
}
