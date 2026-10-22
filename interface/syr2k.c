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
#else
#define ERROR_NAME "SSYR2K"
#endif
#else
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#else
#define ERROR_NAME "CSYR2K"
#endif
#else
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#elif defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#else
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

static int (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#else
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

  blas_arg_t मगा_कातविरान;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
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

  blasint info;
  int uplo;
  int trans;
  int nrowa;

  इलिम्बाग_DEBUG_NAME;

  मगा_कातविरान.n = *N;
  मगा_कातविरान.k = *K;

  मगा_कातविरान.a = (void *)a;
  मगा_कातविरान.b = (void *)b;
  मगा_कातविरान.c = (void *)c;

  मगा_कातविरान.lda = *ldA;
  मगा_कातविरान.ldb = *ldB;
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

  if (मगा_कातविरान.ldc < MAX(1,मगा_कातविरान.n)) info = 12;
  if (मगा_कातविरान.ldb < MAX(1,nrowa))  info =  9;
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
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#else
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
#endif
#if !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#else
	   void *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#else
	   void *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#if !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t मगा_कातविरान;
  int uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
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

  इलिम्बाग_DEBUG_CNAME;

  मगा_कातविरान.n = n;
  मगा_कातविरान.k = k;

  मगा_कातविरान.a = (void *)a;
  मगा_कातविरान.b = (void *)b;
  मगा_कातविरान.c = (void *)c;

  मगा_कातविरान.lda = lda;
  मगा_कातविरान.ldb = ldb;
  मगा_कातविरान.ldc = ldc;

#ifndef COMPLEX
  मगा_कातविरान.alpha = (void *)&alpha;
#else
  मगा_कातविरान.alpha = (void *)alpha;
#endif

#if !defined(COMPLEX) || defined(HEMM)
  मगा_कातविरान.beta  = (void *)&beta;
#else
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

    if (मगा_कातविरान.ldc < MAX(1,मगा_कातविरान.n)) info = 12;
    if (मगा_कातविरान.ldb < MAX(1,nrowa))  info =  9;
    if (मगा_कातविरान.lda < MAX(1,nrowa))  info =  7;
    if (मगा_कातविरान.k < 0)               info =  4;
    if (मगा_कातविरान.n < 0)               info =  3;
    if (trans < 0)                info =  2;
    if (uplo  < 0)                info =  1;
  }

  if (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    मगा_कातविरान.alpha = (void *)CAlpha;
#endif

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

    if (मगा_कातविरान.ldc < MAX(1,मगा_कातविरान.n)) info = 12;
    if (मगा_कातविरान.ldb < MAX(1,nrowa))  info =  9;
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
#if defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#if defined(DYNAMIC_ARCH)
if (strcmp(gotoblas_corename(), "armv9sme") == 0
#if defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   if (order == CblasRowMajor && n == ldc) {
     if (Trans == CblasNoTrans && k == lda && k == ldb) {
      if (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        return;
      }else if (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        return;
      }
     }
     else if (Trans == CblasTrans && n == lda && n ==ldb) {
      if (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        return;
      }else if (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        return;
      }
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
  if (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } else {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  मगा_कातविरान.common = NULL;
  if (मगा_कातविरान.n*मगा_कातविरान.k <1000)
	  मगा_कातविरान.nthreads =1 ;
  else
  मगा_कातविरान.nthreads = नुमेरो_cpu_avail(3);

  if (मगा_कातविरान.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } else {

    syrk_thread(mode, &मगा_कातविरान, NULL, NULL, syr2k[(uplo << 1) | trans ], sa, sb, मगा_कातविरान.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * मगा_कातविरान.n * मगा_कातविरान.k + मगा_कातविरान.n * मगा_कातविरान.n, 2 * मगा_कातविरान.n * मगा_कातविरान.n * मगा_कातविरान.k);

  IDEBUG_END;

  return;
}
