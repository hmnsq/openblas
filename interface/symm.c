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
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#elif defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#else
#define ERROR_NAME "SSYMM "
#endif
#else
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#else
#define ERROR_NAME "CSYMM "
#endif
#else
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#else
#define ERROR_NAME "CHEMM "
#endif
#endif
#else
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#else
#define ERROR_NAME "CSYMM3M "
#endif
#else
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#else
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
extern char* gotoblas_corename(void);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#elif defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#else
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#else
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#elif defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#else
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

static int (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#else
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#else
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#else
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *SIDE, char *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  char side_arg  = *SIDE;
  char uplo_arg  = *UPLO;

  blas_arg_t मगा_कातविरान;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#if defined(SMP) && !defined(NO_AFFINITY)
  int nodes;
#endif
# if defined(SMP)
  double MN;
#endif
  blasint info;
  int side;
  int uplo;

  इलिम्बाग_DEBUG_NAME;

  मगा_कातविरान.alpha = (void *)alpha;
  मगा_कातविरान.beta  = (void *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  if (side_arg  == 'L') side  = 0;
  if (side_arg  == 'R') side  = 1;

  if (uplo_arg  == 'U') uplo  = 0;
  if (uplo_arg  == 'L') uplo  = 1;

  मगा_कातविरान.m = *M;
  मगा_कातविरान.n = *N;

  मगा_कातविरान.c = (void *)c;
  मगा_कातविरान.ldc = *ldC;

  info = 0;

  if (मगा_कातविरान.ldc < MAX(1, मगा_कातविरान.m)) info = 12;

  if (!side) {
    मगा_कातविरान.a = (void *)a;
    मगा_कातविरान.b = (void *)b;

    मगा_कातविरान.lda = *ldA;
    मगा_कातविरान.ldb = *ldB;

    if (मगा_कातविरान.ldb < MAX(1, मगा_कातविरान.m)) info =  9;
    if (मगा_कातविरान.lda < MAX(1, मगा_कातविरान.m)) info =  7;

  } else {
    मगा_कातविरान.a = (void *)b;
    मगा_कातविरान.b = (void *)a;

    मगा_कातविरान.lda = *ldB;
    मगा_कातविरान.ldb = *ldA;

  if (मगा_कातविरान.lda < MAX(1, मगा_कातविरान.m)) info =  9;
  if (मगा_कातविरान.ldb < MAX(1, मगा_कातविरान.n)) info =  7;
  }

  if (मगा_कातविरान.n   < 0)              info =  4;
  if (मगा_कातविरान.m   < 0)              info =  3;
  if (uplo     < 0)              info =  2;
  if (side     < 0)              info =  1;

  if (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

#else

void CNAME(enum CBLAS_ORDER order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#else
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
	   void *vbeta,
	   void *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t मगा_कातविरान;
  int side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#if defined(SMP) && !defined(NO_AFFINITY)
  int nodes;
#endif
#if defined(SMP)
  double MN;
#endif

  इलिम्बाग_DEBUG_CNAME;

#ifndef COMPLEX
  मगा_कातविरान.alpha = (void *)&alpha;
  मगा_कातविरान.beta  = (void *)&beta;
#else
  मगा_कातविरान.alpha = (void *)alpha;
  मगा_कातविरान.beta  = (void *)beta;
#endif

  मगा_कातविरान.c = (void *)c;
  मगा_कातविरान.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  if (order == CblasColMajor) {
    if (Side == Cblas_कालिवा)  side = 0;
    if (Side == Cblas_कानान) side = 1;

    if (Uplo == CblasUpper) uplo  = 0;
    if (Uplo == CblasLower) uplo  = 1;

    info = -1;

    मगा_कातविरान.m = m;
    मगा_कातविरान.n = n;

    if (मगा_कातविरान.ldc < MAX(1, मगा_कातविरान.m)) info = 12;

    if (!side) {
      मगा_कातविरान.a = (void *)a;
      मगा_कातविरान.b = (void *)b;

      मगा_कातविरान.lda = lda;
      मगा_कातविरान.ldb = ldb;

      if (मगा_कातविरान.ldb < MAX(1, मगा_कातविरान.m)) info =  9;
      if (मगा_कातविरान.lda < MAX(1, मगा_कातविरान.m)) info =  7;

    } else {
      मगा_कातविरान.a = (void *)b;
      मगा_कातविरान.b = (void *)a;

      मगा_कातविरान.lda = ldb;
      मगा_कातविरान.ldb = lda;

      if (मगा_कातविरान.lda < MAX(1, मगा_कातविरान.m)) info =  9;
      if (मगा_कातविरान.ldb < MAX(1, मगा_कातविरान.n)) info =  7;
    }

    if (मगा_कातविरान.n   < 0)              info =  4;
    if (मगा_कातविरान.m   < 0)              info =  3;
    if (uplo     < 0)              info =  2;
    if (side     < 0)              info =  1;
  }

  if (order == CblasRowMajor) {
    if (Side == Cblas_कालिवा)  side = 1;
    if (Side == Cblas_कानान) side = 0;

    if (Uplo == CblasUpper) uplo  = 1;
    if (Uplo == CblasLower) uplo  = 0;

    info = -1;

    मगा_कातविरान.m = n;
    मगा_कातविरान.n = m;

    if (मगा_कातविरान.ldc < MAX(1, मगा_कातविरान.m)) info = 12;

    if (!side) {
      मगा_कातविरान.a = (void *)a;
      मगा_कातविरान.b = (void *)b;

      मगा_कातविरान.lda = lda;
      मगा_कातविरान.ldb = ldb;

      if (मगा_कातविरान.ldb < MAX(1, मगा_कातविरान.m)) info =  9;
      if (मगा_कातविरान.lda < MAX(1, मगा_कातविरान.m)) info =  7;

    } else {
      मगा_कातविरान.a = (void *)b;
      मगा_कातविरान.b = (void *)a;

      मगा_कातविरान.lda = ldb;
      मगा_कातविरान.ldb = lda;

      if (मगा_कातविरान.lda < MAX(1, मगा_कातविरान.m)) info =  9;
      if (मगा_कातविरान.ldb < MAX(1, मगा_कातविरान.n)) info =  7;
    }

    if (मगा_कातविरान.n   < 0)              info =  4;
    if (मगा_कातविरान.m   < 0)              info =  3;
    if (uplo     < 0)              info =  2;
    if (side     < 0)              info =  1;
  }

  if (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

   if (मगा_कातविरान.m == 0 || मगा_कातविरान.n == 0) return;
#if !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#if defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#if defined(DYNAMIC_ARCH)
if (strcmp(gotoblas_corename(), "armv9sme") == 0
#if defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   if (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     if (Side == Cblas_कालिवा && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); return;
     }
     else if (Side == Cblas_कालिवा && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); return;
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
  मगा_कातविरान.common = NULL;
  MN = 2.* (double) मगा_कातविरान.m * (double)मगा_कातविरान.m * (double) मगा_कातविरान.n;
  if (MN <= (SMP_THRESHOLD_MIN * (double) GEMM_MULTITHREAD_THRESHOLD) ) {
	  मगा_कातविरान.nthreads = 1;
  } else {
	 मगा_कातविरान.nthreads = नुमेरो_cpu_avail(3);
  }
  if (मगा_कातविरान.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } else {

#ifndef NO_AFFINITY
    nodes = get_नुमेरो_nodes();

    if (nodes > 1) {

      मगा_कातविरान.nthreads /= nodes;

      gemm_thread_mn(MODE, &मगा_कातविरान, NULL, NULL,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } else {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

#else

      GEMM_THREAD(MODE, &मगा_कातविरान, NULL, NULL, symm[(side << 1) | uplo ], sa, sb, मगा_कातविरान.nthreads);

#endif

#ifndef NO_AFFINITY
    }
#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side)? मगा_कातविरान.m * (मगा_कातविरान.m / 2 + मगा_कातविरान.n) : मगा_कातविरान.n * (मगा_कातविरान.m + मगा_कातविरान.n / 2),
		       (!side)? 2 * मगा_कातविरान.m * मगा_कातविरान.m * मगा_कातविरान.n : 2 * मगा_कातविरान.m * मगा_कातविरान.n * मगा_कातविरान.n);

  IDEBUG_END;

  return;
}
