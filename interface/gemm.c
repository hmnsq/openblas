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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <stdlib.h>
#അന്തര്ഭാവ <stdbool.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#അല്ലെങ്കിലെങ്കിൽ defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#അല്ലെങ്കിൽ
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#അല്ലെങ്കിലെങ്കിൽ defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#അല്ലെങ്കിൽ
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#അല്ലെങ്കിൽ
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#എങ്കിൽ !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#എങ്കിൽ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#അല്ലെങ്കിൽ
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#എങ്കിൽ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#എങ്കിൽ defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#അല്ലെങ്കിൽ
#define USE_SMALL_MATRIX_OPT 0
#endif

#എങ്കിൽ USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((ശൂന്യം *)(table[idx]))
#അല്ലെങ്കിൽ
#define SMALL_KERNEL_ADDR(table, idx) ((ശൂന്യം *)(*(uintptr_t *)((അക്ഷരം *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
സ്ഥിതം size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


സ്ഥിതം size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#അല്ലെങ്കിൽ

സ്ഥിതം size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

സ്ഥിതം size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#എങ്കിൽ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
സ്ഥിതം പൂർണ്ണസംഖ്യ openblas_amxtile_permission = 0;
സ്ഥിതം പൂർണ്ണസംഖ്യ init_amxtile_permission() {
  long status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  എങ്കിൽ (status != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    തിരികെ -1;
  }
  openblas_amxtile_permission = 1;
  തിരികെ 0;
}
#endif

#ifdef DYNAMIC_ARCH
ബാഹ്യ അക്ഷരം* gotoblas_corename(ശൂന്യം);
#endif

#ifdef SMP
#എങ്കിൽ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
സ്ഥിതം വരിയിൽ പൂർണ്ണസംഖ്യ get_gemm_optimal_nthreads_neoversev1(ഇരട്ട MNK, പൂർണ്ണസംഖ്യ ncpu) {
  തിരികെ
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

#എങ്കിൽ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
സ്ഥിതം വരിയിൽ പൂർണ്ണസംഖ്യ get_gemm_optimal_nthreads_neoversev2(ഇരട്ട MNK, പൂർണ്ണസംഖ്യ ncpu) {
  തിരികെ
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

സ്ഥിതം വരിയിൽ പൂർണ്ണസംഖ്യ get_gemm_optimal_nthreads(ഇരട്ട MNK) {
  പൂർണ്ണസംഖ്യ ncpu = അക്കം_cpu_avail(3);
#എങ്കിൽ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  തിരികെ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#അല്ലെങ്കിലെങ്കിൽ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  തിരികെ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#അല്ലെങ്കിലെങ്കിൽ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  എങ്കിൽ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    തിരികെ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  എങ്കിൽ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    തിരികെ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  എങ്കിൽ ( MNK <= (SMP_THRESHOLD_MIN  * (ഇരട്ട) GEMM_MULTITHREAD_THRESHOLD) ) {
    തിരികെ 1;
  }
  അല്ലെങ്കിൽ {
    എങ്കിൽ (MNK/ncpu < SMP_THRESHOLD_MIN*(ഇരട്ട)GEMM_MULTITHREAD_THRESHOLD) {
      തിരികെ MNK/(SMP_THRESHOLD_MIN*(ഇരട്ട)GEMM_MULTITHREAD_THRESHOLD);
    }
    അല്ലെങ്കിൽ {
      തിരികെ ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *TRANSA, അക്ഷരം *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  പൂർണ്ണസംഖ്യ transa, transb, nrowa, nrowb;
  blasint info;

  അക്ഷരം transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  ഇരട്ട MNK;
#എങ്കിൽ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#എങ്കിൽ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  പൂർണ്ണസംഖ്യ nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (ശൂന്യം *)a;
  args.b = (ശൂന്യം *)b;
  args.c = (ശൂന്യം *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (ശൂന്യം *)alpha;
  args.beta  = (ശൂന്യം *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  എങ്കിൽ (transA == 'N') transa = 0;
  എങ്കിൽ (transA == 'T') transa = 1;
#ifndef COMPLEX
  എങ്കിൽ (transA == 'R') transa = 0;
  എങ്കിൽ (transA == 'C') transa = 1;
#അല്ലെങ്കിൽ
  എങ്കിൽ (transA == 'R') transa = 2;
  എങ്കിൽ (transA == 'C') transa = 3;
#endif

  എങ്കിൽ (transB == 'N') transb = 0;
  എങ്കിൽ (transB == 'T') transb = 1;
#ifndef COMPLEX
  എങ്കിൽ (transB == 'R') transb = 0;
  എങ്കിൽ (transB == 'C') transb = 1;
#അല്ലെങ്കിൽ
  എങ്കിൽ (transB == 'R') transb = 2;
  എങ്കിൽ (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  എങ്കിൽ (transa & 1) nrowa = args.k;
  nrowb = args.k;
  എങ്കിൽ (transb & 1) nrowb = args.n;

  info = 0;

  എങ്കിൽ (args.ldc < args.m) info = 13;
  എങ്കിൽ (args.ldb < nrowb)  info = 10;
  എങ്കിൽ (args.lda < nrowa)  info =  8;
  എങ്കിൽ (args.k < 0)        info =  5;
  എങ്കിൽ (args.n < 0)        info =  4;
  എങ്കിൽ (args.m < 0)        info =  3;
  എങ്കിൽ (transb < 0)        info =  2;
  എങ്കിൽ (transa < 0)        info =  1;

  എങ്കിൽ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_TRANSPOSE TransA, എണ്ണൽ CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#അല്ലെങ്കിൽ
	   ശൂന്യം *valpha,
	   ശൂന്യം *va, blasint lda,
	   ശൂന്യം *vb, blasint ldb,
	   ശൂന്യം *vbeta,
	   ശൂന്യം *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  പൂർണ്ണസംഖ്യ transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  ഇരട്ട MNK;
#എങ്കിൽ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#എങ്കിൽ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  പൂർണ്ണസംഖ്യ nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (ശൂന്യം *)&alpha;
  args.beta  = (ശൂന്യം *)&beta;
#അല്ലെങ്കിൽ
  args.alpha = (ശൂന്യം *)alpha;
  args.beta  = (ശൂന്യം *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  എങ്കിൽ (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (ശൂന്യം *)a;
    args.b = (ശൂന്യം *)b;
    args.c = (ശൂന്യം *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    എങ്കിൽ (TransA == CblasNoTrans)     transa = 0;
    എങ്കിൽ (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    എങ്കിൽ (TransA == CblasConjNoTrans) transa = 0;
    എങ്കിൽ (TransA == CblasConjTrans)   transa = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (TransA == CblasConjNoTrans) transa = 2;
    എങ്കിൽ (TransA == CblasConjTrans)   transa = 3;
#endif
    എങ്കിൽ (TransB == CblasNoTrans)     transb = 0;
    എങ്കിൽ (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    എങ്കിൽ (TransB == CblasConjNoTrans) transb = 0;
    എങ്കിൽ (TransB == CblasConjTrans)   transb = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (TransB == CblasConjNoTrans) transb = 2;
    എങ്കിൽ (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    എങ്കിൽ (transa & 1) nrowa = args.k;
    nrowb = args.k;
    എങ്കിൽ (transb & 1) nrowb = args.n;

    info = -1;

    എങ്കിൽ (args.ldc < args.m) info = 13;
    എങ്കിൽ (args.ldb < nrowb)  info = 10;
    എങ്കിൽ (args.lda < nrowa)  info =  8;
    എങ്കിൽ (args.k < 0)        info =  5;
    എങ്കിൽ (args.n < 0)        info =  4;
    എങ്കിൽ (args.m < 0)        info =  3;
    എങ്കിൽ (transb < 0)        info =  2;
    എങ്കിൽ (transa < 0)        info =  1;
  }

  എങ്കിൽ (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (ശൂന്യം *)b;
    args.b = (ശൂന്യം *)a;
    args.c = (ശൂന്യം *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    എങ്കിൽ (TransB == CblasNoTrans)     transa = 0;
    എങ്കിൽ (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    എങ്കിൽ (TransB == CblasConjNoTrans) transa = 0;
    എങ്കിൽ (TransB == CblasConjTrans)   transa = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (TransB == CblasConjNoTrans) transa = 2;
    എങ്കിൽ (TransB == CblasConjTrans)   transa = 3;
#endif
    എങ്കിൽ (TransA == CblasNoTrans)     transb = 0;
    എങ്കിൽ (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    എങ്കിൽ (TransA == CblasConjNoTrans) transb = 0;
    എങ്കിൽ (TransA == CblasConjTrans)   transb = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (TransA == CblasConjNoTrans) transb = 2;
    എങ്കിൽ (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    എങ്കിൽ (transa & 1) nrowa = args.k;
    nrowb = args.k;
    എങ്കിൽ (transb & 1) nrowb = args.n;

    info = -1;

    എങ്കിൽ (args.ldc < args.m) info = 13;
    എങ്കിൽ (args.ldb < nrowb)  info = 10;
    എങ്കിൽ (args.lda < nrowa)  info =  8;
    എങ്കിൽ (args.k < 0)        info =  5;
    എങ്കിൽ (args.n < 0)        info =  4;
    എങ്കിൽ (args.m < 0)        info =  3;
    എങ്കിൽ (transb < 0)        info =  2;
    എങ്കിൽ (transa < 0)        info =  1;

  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }


  എങ്കിൽ ((args.m == 0) || (args.n == 0)) തിരികെ;
#എങ്കിൽ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#എങ്കിൽ defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#എങ്കിൽ defined(DYNAMIC_ARCH)
  എങ്കിൽ (support_avx512() )
#endif
  എങ്കിൽ (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        തിരികെ;
  }
#endif
#എങ്കിൽ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#എങ്കിൽ defined(DYNAMIC_ARCH)
എങ്കിൽ (strcmp(gotoblas_corename(), "armv9sme") == 0
#എങ്കിൽ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  എങ്കിൽ (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        തിരികെ;
  }
അല്ലെങ്കിൽ എങ്കിൽ (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        തിരികെ;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#എങ്കിൽ !defined(BFLOAT16)  && !defined(HFLOAT16)
#എങ്കിൽ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#എങ്കിൽ defined(DYNAMIC_ARCH)
എങ്കിൽ (strcmp(gotoblas_corename(), "armv9sme") == 0
#എങ്കിൽ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
അക്ഷരം* TA,*TB;
      എങ്കിൽ (transa & 1) 
        TA = "T";
	അല്ലെങ്കിൽ 
	TA= "N";
      എങ്കിൽ (transb & 1) 
        TB = "T";
	അല്ലെങ്കിൽ 
	TB= "N";
#ifndef COMPLEX
      എങ്കിൽ (transa == 3) 
	TA= "T";
      എങ്കിൽ (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#അല്ലെങ്കിൽ
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#അല്ലെങ്കിൽ
      എങ്കിൽ (transa == 2) 
	TA= "R";
      എങ്കിൽ (transb == 2) 
	TB= "R";
      എങ്കിൽ (transa == 3) 
	TA= "C";
      എങ്കിൽ (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
ദശാംശം _Complex c_al={al[0],al[1]};
ദശാംശം _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#അല്ലെങ്കിൽ
ഇരട്ട _Complex c_al={al[0],al[1]};
ഇരട്ട _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	തിരികെ;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#എങ്കിൽ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#എങ്കിൽ defined(DYNAMIC_ARCH)
  എങ്കിൽ (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    തിരികെ;
  }
#endif
#എങ്കിൽ !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  എങ്കിൽ (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    തിരികെ;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  എങ്കിൽ ((args.m == 0) || (args.n == 0)) തിരികെ;

#എങ്കിൽ 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#എങ്കിൽ (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#അല്ലെങ്കിൽ
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#എങ്കിൽ (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#അല്ലെങ്കിൽ
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#എങ്കിൽ defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#എങ്കിൽ defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#അല്ലെങ്കിൽ
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  സ്ഥിതം സ്ഥിരം അക്ഷരം gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  എങ്കിൽ (args.k != 0) {
    എങ്കിൽ (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      പൂർണ്ണസംഖ്യ trans_idx = transa | ((transb & 2) << 1);
      അക്ഷരം NT = gemv_trans[trans_idx];
      എങ്കിൽ (transa & 1) {
        m = args.k;
        n = args.m;
      }
      എങ്കിൽ (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      എങ്കിൽ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        തിരികെ;
      }
    }
    എങ്കിൽ (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      പൂർണ്ണസംഖ്യ trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      അക്ഷരം NT = gemv_trans[trans_idx];
      എങ്കിൽ (transa & 1) {
        inc_x = 1;
      }
      എങ്കിൽ (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      എങ്കിൽ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        തിരികെ;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#എങ്കിൽ USE_SMALL_MATRIX_OPT
എങ്കിൽ (*(FLOAT*)args.alpha != 0.0) {
#എങ്കിൽ !defined(COMPLEX)
  എങ്കിൽ(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  എങ്കിൽ(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }അല്ലെങ്കിൽ{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  തിരികെ;
  }
#അല്ലെങ്കിൽ
  എങ്കിൽ(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  എങ്കിൽ(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }അല്ലെങ്കിൽ{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  തിരികെ;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#എങ്കിൽ defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#അല്ലെങ്കിൽ
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#എങ്കിൽ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (ഇരട്ട) args.m * (ഇരട്ട) args.n * (ഇരട്ട) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = ശൂന്യ;

 എങ്കിൽ (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#ifdef SMP

  } അല്ലെങ്കിൽ {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_അക്കം_nodes();

      എങ്കിൽ ((nodes > 1) && get_കണ്ണി_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, ശൂന്യ, ശൂന്യ, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } അല്ലെങ്കിൽ {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#അല്ലെങ്കിൽ

	GEMM_THREAD(mode, &args, ശൂന്യ, ശൂന്യ, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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

  തിരികെ;
}
