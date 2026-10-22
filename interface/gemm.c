/*********************************************************************/
/* Copyright 2024, 2025 The OpenBLAS Project                         */
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
#include <stdlib.h>
#include <stdbool.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#elif defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#không_thì
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#elif defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#không_thì
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#không_thì
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#không_thì
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#không_thì
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#không_thì
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

static số_nguyên (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#nếu !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#nếu defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#không_thì
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#nếu defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#nếu defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#không_thì
#define USE_SMALL_MATRIX_OPT 0
#endif

#nếu USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(table[idx]))
#không_thì
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(*(uintptr_t *)((char *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
static size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


static size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#không_thì

static size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

static size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#nếu defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
static số_nguyên openblas_amxtile_permission = 0;
static số_nguyên init_amxtile_permission() {
  long status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  nếu (status != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    trả_về -1;
  }
  openblas_amxtile_permission = 1;
  trả_về 0;
}
#endif

#ifdef DYNAMIC_ARCH
extern char* gotoblas_corename(void);
#endif

#ifdef SMP
#nếu defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
static inline số_nguyên get_gemm_optimal_nthreads_neoversev1(double MNK, số_nguyên ncpu) {
  trả_về
      MNK < 262144L    ? 1
    : MNK < 1124864L   ? MIN(ncpu, 6)
    : MNK < 7880599L   ? MIN(ncpu, 12)
    : MNK < 17173512L  ? MIN(ncpu, 16)
    : MNK < 33386248L  ? MIN(ncpu, 20)
    : MNK < 57066625L  ? MIN(ncpu, 24)
    : MNK < 91733851L  ? MIN(ncpu, 32)
    : MNK < 265847707L ? MIN(ncpu, 40)
    : MNK < 458314011L ? MIN(ncpu, 48)
    : MNK < 729000000L ? MIN(ncpu, 56)
    : ncpu;
}
#endif

#nếu defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
static inline số_nguyên get_gemm_optimal_nthreads_neoversev2(double MNK, số_nguyên ncpu) {
  trả_về
      MNK < 125000L     ? 1
    : MNK < 1092727L    ? MIN(ncpu, 6)
    : MNK < 2628072L    ? MIN(ncpu, 8)
    : MNK < 8000000L    ? MIN(ncpu, 12)
    : MNK < 20346417L   ? MIN(ncpu, 16)
    : MNK < 57066625L   ? MIN(ncpu, 24)
    : MNK < 91125000L   ? MIN(ncpu, 28)
    : MNK < 238328000L  ? MIN(ncpu, 40)
    : MNK < 454756609L  ? MIN(ncpu, 48)
    : MNK < 857375000L  ? MIN(ncpu, 56)
    : MNK < 1073741824L ? MIN(ncpu, 64)
    : ncpu;
}
#endif

static inline số_nguyên get_gemm_optimal_nthreads(double MNK) {
  số_nguyên ncpu = num_cpu_avail(3);
#nếu defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  trả_về get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#elif defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  trả_về get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#elif defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  nếu (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    trả_về get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  nếu (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    trả_về get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  nếu ( MNK <= (SMP_THRESHOLD_MIN  * (double) GEMM_MULTITHREAD_THRESHOLD) ) {
    trả_về 1;
  }
  không_thì {
    nếu (MNK/ncpu < SMP_THRESHOLD_MIN*(double)GEMM_MULTITHREAD_THRESHOLD) {
      trả_về MNK/(SMP_THRESHOLD_MIN*(double)GEMM_MULTITHREAD_THRESHOLD);
    }
    không_thì {
      trả_về ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *TRANSA, char *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  số_nguyên transa, transb, nrowa, nrowb;
  blasint info;

  char transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  double MNK;
#nếu defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  số_nguyên mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  số_nguyên mode  =  BLAS_DOUBLE  | BLAS_REAL;
#không_thì
  số_nguyên mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#không_thì
#ifdef XDOUBLE
  số_nguyên mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  số_nguyên mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#không_thì
  số_nguyên mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#nếu defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  số_nguyên nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
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
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  nếu (transA == 'N') transa = 0;
  nếu (transA == 'T') transa = 1;
#ifndef COMPLEX
  nếu (transA == 'R') transa = 0;
  nếu (transA == 'C') transa = 1;
#không_thì
  nếu (transA == 'R') transa = 2;
  nếu (transA == 'C') transa = 3;
#endif

  nếu (transB == 'N') transb = 0;
  nếu (transB == 'T') transb = 1;
#ifndef COMPLEX
  nếu (transB == 'R') transb = 0;
  nếu (transB == 'C') transb = 1;
#không_thì
  nếu (transB == 'R') transb = 2;
  nếu (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  nếu (transa & 1) nrowa = args.k;
  nrowb = args.k;
  nếu (transb & 1) nrowb = args.n;

  info = 0;

  nếu (args.ldc < args.m) info = 13;
  nếu (args.ldb < nrowb)  info = 10;
  nếu (args.lda < nrowa)  info =  8;
  nếu (args.k < 0)        info =  5;
  nếu (args.n < 0)        info =  4;
  nếu (args.m < 0)        info =  3;
  nếu (transb < 0)        info =  2;
  nếu (transa < 0)        info =  1;

  nếu (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#không_thì

void CNAME(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#không_thì
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

  blas_arg_t args;
  số_nguyên transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  double MNK;
#nếu defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  số_nguyên mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  số_nguyên mode  =  BLAS_DOUBLE  | BLAS_REAL;
#không_thì
  số_nguyên mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#không_thì
#ifdef XDOUBLE
  số_nguyên mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  số_nguyên mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#không_thì
  số_nguyên mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#nếu defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  số_nguyên nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (void *)&alpha;
  args.beta  = (void *)&beta;
#không_thì
  args.alpha = (void *)alpha;
  args.beta  = (void *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  nếu (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (void *)a;
    args.b = (void *)b;
    args.c = (void *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    nếu (TransA == CblasNoTrans)     transa = 0;
    nếu (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    nếu (TransA == CblasConjNoTrans) transa = 0;
    nếu (TransA == CblasConjTrans)   transa = 1;
#không_thì
    nếu (TransA == CblasConjNoTrans) transa = 2;
    nếu (TransA == CblasConjTrans)   transa = 3;
#endif
    nếu (TransB == CblasNoTrans)     transb = 0;
    nếu (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    nếu (TransB == CblasConjNoTrans) transb = 0;
    nếu (TransB == CblasConjTrans)   transb = 1;
#không_thì
    nếu (TransB == CblasConjNoTrans) transb = 2;
    nếu (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    nếu (transa & 1) nrowa = args.k;
    nrowb = args.k;
    nếu (transb & 1) nrowb = args.n;

    info = -1;

    nếu (args.ldc < args.m) info = 13;
    nếu (args.ldb < nrowb)  info = 10;
    nếu (args.lda < nrowa)  info =  8;
    nếu (args.k < 0)        info =  5;
    nếu (args.n < 0)        info =  4;
    nếu (args.m < 0)        info =  3;
    nếu (transb < 0)        info =  2;
    nếu (transa < 0)        info =  1;
  }

  nếu (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (void *)b;
    args.b = (void *)a;
    args.c = (void *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    nếu (TransB == CblasNoTrans)     transa = 0;
    nếu (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    nếu (TransB == CblasConjNoTrans) transa = 0;
    nếu (TransB == CblasConjTrans)   transa = 1;
#không_thì
    nếu (TransB == CblasConjNoTrans) transa = 2;
    nếu (TransB == CblasConjTrans)   transa = 3;
#endif
    nếu (TransA == CblasNoTrans)     transb = 0;
    nếu (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    nếu (TransA == CblasConjNoTrans) transb = 0;
    nếu (TransA == CblasConjTrans)   transb = 1;
#không_thì
    nếu (TransA == CblasConjNoTrans) transb = 2;
    nếu (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    nếu (transa & 1) nrowa = args.k;
    nrowb = args.k;
    nếu (transb & 1) nrowb = args.n;

    info = -1;

    nếu (args.ldc < args.m) info = 13;
    nếu (args.ldb < nrowb)  info = 10;
    nếu (args.lda < nrowa)  info =  8;
    nếu (args.k < 0)        info =  5;
    nếu (args.n < 0)        info =  4;
    nếu (args.m < 0)        info =  3;
    nếu (transb < 0)        info =  2;
    nếu (transa < 0)        info =  1;

  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }


  nếu ((args.m == 0) || (args.n == 0)) trả_về;
#nếu !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#nếu defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#nếu defined(DYNAMIC_ARCH)
  nếu (support_avx512() )
#endif
  nếu (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        trả_về;
  }
#endif
#nếu defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#nếu defined(DYNAMIC_ARCH)
nếu (strcmp(gotoblas_corename(), "armv9sme") == 0
#nếu defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  nếu (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        trả_về;
  }
không_thì nếu (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        trả_về;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#nếu !defined(BFLOAT16)  && !defined(HFLOAT16)
#nếu defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#nếu defined(DYNAMIC_ARCH)
nếu (strcmp(gotoblas_corename(), "armv9sme") == 0
#nếu defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
char* TA,*TB;
      nếu (transa & 1) 
        TA = "T";
	không_thì 
	TA= "N";
      nếu (transb & 1) 
        TB = "T";
	không_thì 
	TB= "N";
#ifndef COMPLEX
      nếu (transa == 3) 
	TA= "T";
      nếu (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#không_thì
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#không_thì
      nếu (transa == 2) 
	TA= "R";
      nếu (transb == 2) 
	TB= "R";
      nếu (transa == 3) 
	TA= "C";
      nếu (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
float _Complex c_al={al[0],al[1]};
float _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#không_thì
double _Complex c_al={al[0],al[1]};
double _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	trả_về;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#nếu defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#nếu defined(DYNAMIC_ARCH)
  nếu (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    trả_về;
  }
#endif
#nếu !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  nếu (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    trả_về;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  nếu ((args.m == 0) || (args.n == 0)) trả_về;

#nếu 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#nếu (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#không_thì
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#nếu (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#không_thì
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#nếu defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#nếu defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#không_thì
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  static const char gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  nếu (args.k != 0) {
    nếu (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      số_nguyên trans_idx = transa | ((transb & 2) << 1);
      char NT = gemv_trans[trans_idx];
      nếu (transa & 1) {
        m = args.k;
        n = args.m;
      }
      nếu (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      nếu (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        trả_về;
      }
    }
    nếu (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      số_nguyên trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      char NT = gemv_trans[trans_idx];
      nếu (transa & 1) {
        inc_x = 1;
      }
      nếu (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      nếu (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        trả_về;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#nếu USE_SMALL_MATRIX_OPT
nếu (*(FLOAT*)args.alpha != 0.0) {
#nếu !defined(COMPLEX)
  nếu(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  nếu(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }không_thì{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  trả_về;
  }
#không_thì
  nếu(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  nếu(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }không_thì{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  trả_về;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#nếu defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#không_thì
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#nếu defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (double) args.m * (double) args.n * (double) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = NULL;

 nếu (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } không_thì {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_num_nodes();

      nếu ((nodes > 1) && get_node_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, NULL, NULL, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } không_thì {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, NULL, NULL, sa, sb, 0);

#không_thì

	GEMM_THREAD(mode, &args, NULL, NULL, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3
#ifndef NO_AFFINITY
      }
#endif
#endif

#endif

#ifdef SMP
  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n, 2 * args.m * args.n * args.k);

  IDEBUG_END;

  trả_về;
}
