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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <stdlib.h>
#அந்தர்பாவ <stdbool.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#அல்லதெனில் defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#இல்லையெனில்
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#அல்லதெனில் defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#இல்லையெனில்
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#இல்லையெனில்
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#இல்லையெனில்
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#இல்லையெனில்
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

நிலையான முழு (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#எனில் !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#இல்லையெனில்
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#எனில் defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#இல்லையெனில்
#define USE_SMALL_MATRIX_OPT 0
#endif

#எனில் USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((வெறுமை *)(table[idx]))
#இல்லையெனில்
#define SMALL_KERNEL_ADDR(table, idx) ((வெறுமை *)(*(uintptr_t *)((எழுத்து *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
நிலையான size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


நிலையான size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (முழு (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (முழு (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#இல்லையெனில்

நிலையான size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

நிலையான size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (முழு (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (முழு (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#எனில் defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
நிலையான முழு openblas_amxtile_permission = 0;
நிலையான முழு init_amxtile_permission() {
  long status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  எனில் (status != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    திரும்பு -1;
  }
  openblas_amxtile_permission = 1;
  திரும்பு 0;
}
#endif

#ifdef DYNAMIC_ARCH
வெளி எழுத்து* gotoblas_corename(வெறுமை);
#endif

#ifdef SMP
#எனில் defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
நிலையான உள்வரி முழு get_gemm_optimal_nthreads_neoversev1(இரட்டை MNK, முழு ncpu) {
  திரும்பு
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

#எனில் defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
நிலையான உள்வரி முழு get_gemm_optimal_nthreads_neoversev2(இரட்டை MNK, முழு ncpu) {
  திரும்பு
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

நிலையான உள்வரி முழு get_gemm_optimal_nthreads(இரட்டை MNK) {
  முழு ncpu = இலக்கம்_cpu_avail(3);
#எனில் defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  திரும்பு get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#அல்லதெனில் defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  திரும்பு get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#அல்லதெனில் defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  எனில் (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    திரும்பு get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  எனில் (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    திரும்பு get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  எனில் ( MNK <= (SMP_THRESHOLD_MIN  * (இரட்டை) GEMM_MULTITHREAD_THRESHOLD) ) {
    திரும்பு 1;
  }
  இல்லையெனில் {
    எனில் (MNK/ncpu < SMP_THRESHOLD_MIN*(இரட்டை)GEMM_MULTITHREAD_THRESHOLD) {
      திரும்பு MNK/(SMP_THRESHOLD_MIN*(இரட்டை)GEMM_MULTITHREAD_THRESHOLD);
    }
    இல்லையெனில் {
      திரும்பு ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *TRANSA, எழுத்து *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  முழு transa, transb, nrowa, nrowb;
  blasint info;

  எழுத்து transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  இரட்டை MNK;
#எனில் defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#எனில் defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  முழு nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (வெறுமை *)a;
  args.b = (வெறுமை *)b;
  args.c = (வெறுமை *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  எனில் (transA == 'N') transa = 0;
  எனில் (transA == 'T') transa = 1;
#ifndef COMPLEX
  எனில் (transA == 'R') transa = 0;
  எனில் (transA == 'C') transa = 1;
#இல்லையெனில்
  எனில் (transA == 'R') transa = 2;
  எனில் (transA == 'C') transa = 3;
#endif

  எனில் (transB == 'N') transb = 0;
  எனில் (transB == 'T') transb = 1;
#ifndef COMPLEX
  எனில் (transB == 'R') transb = 0;
  எனில் (transB == 'C') transb = 1;
#இல்லையெனில்
  எனில் (transB == 'R') transb = 2;
  எனில் (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  எனில் (transa & 1) nrowa = args.k;
  nrowb = args.k;
  எனில் (transb & 1) nrowb = args.n;

  info = 0;

  எனில் (args.ldc < args.m) info = 13;
  எனில் (args.ldb < nrowb)  info = 10;
  எனில் (args.lda < nrowa)  info =  8;
  எனில் (args.k < 0)        info =  5;
  எனில் (args.n < 0)        info =  4;
  எனில் (args.m < 0)        info =  3;
  எனில் (transb < 0)        info =  2;
  எனில் (transa < 0)        info =  1;

  எனில் (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_TRANSPOSE TransA, எண்ணிக்கை CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#இல்லையெனில்
	   வெறுமை *valpha,
	   வெறுமை *va, blasint lda,
	   வெறுமை *vb, blasint ldb,
	   வெறுமை *vbeta,
	   வெறுமை *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  முழு transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  இரட்டை MNK;
#எனில் defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#எனில் defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  முழு nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (வெறுமை *)&alpha;
  args.beta  = (வெறுமை *)&beta;
#இல்லையெனில்
  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  எனில் (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (வெறுமை *)a;
    args.b = (வெறுமை *)b;
    args.c = (வெறுமை *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    எனில் (TransA == CblasNoTrans)     transa = 0;
    எனில் (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    எனில் (TransA == CblasConjNoTrans) transa = 0;
    எனில் (TransA == CblasConjTrans)   transa = 1;
#இல்லையெனில்
    எனில் (TransA == CblasConjNoTrans) transa = 2;
    எனில் (TransA == CblasConjTrans)   transa = 3;
#endif
    எனில் (TransB == CblasNoTrans)     transb = 0;
    எனில் (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    எனில் (TransB == CblasConjNoTrans) transb = 0;
    எனில் (TransB == CblasConjTrans)   transb = 1;
#இல்லையெனில்
    எனில் (TransB == CblasConjNoTrans) transb = 2;
    எனில் (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    எனில் (transa & 1) nrowa = args.k;
    nrowb = args.k;
    எனில் (transb & 1) nrowb = args.n;

    info = -1;

    எனில் (args.ldc < args.m) info = 13;
    எனில் (args.ldb < nrowb)  info = 10;
    எனில் (args.lda < nrowa)  info =  8;
    எனில் (args.k < 0)        info =  5;
    எனில் (args.n < 0)        info =  4;
    எனில் (args.m < 0)        info =  3;
    எனில் (transb < 0)        info =  2;
    எனில் (transa < 0)        info =  1;
  }

  எனில் (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (வெறுமை *)b;
    args.b = (வெறுமை *)a;
    args.c = (வெறுமை *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    எனில் (TransB == CblasNoTrans)     transa = 0;
    எனில் (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    எனில் (TransB == CblasConjNoTrans) transa = 0;
    எனில் (TransB == CblasConjTrans)   transa = 1;
#இல்லையெனில்
    எனில் (TransB == CblasConjNoTrans) transa = 2;
    எனில் (TransB == CblasConjTrans)   transa = 3;
#endif
    எனில் (TransA == CblasNoTrans)     transb = 0;
    எனில் (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    எனில் (TransA == CblasConjNoTrans) transb = 0;
    எனில் (TransA == CblasConjTrans)   transb = 1;
#இல்லையெனில்
    எனில் (TransA == CblasConjNoTrans) transb = 2;
    எனில் (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    எனில் (transa & 1) nrowa = args.k;
    nrowb = args.k;
    எனில் (transb & 1) nrowb = args.n;

    info = -1;

    எனில் (args.ldc < args.m) info = 13;
    எனில் (args.ldb < nrowb)  info = 10;
    எனில் (args.lda < nrowa)  info =  8;
    எனில் (args.k < 0)        info =  5;
    எனில் (args.n < 0)        info =  4;
    எனில் (args.m < 0)        info =  3;
    எனில் (transb < 0)        info =  2;
    எனில் (transa < 0)        info =  1;

  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }


  எனில் ((args.m == 0) || (args.n == 0)) திரும்பு;
#எனில் !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#எனில் defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
  எனில் (support_avx512() )
#endif
  எனில் (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        திரும்பு;
  }
#endif
#எனில் defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
எனில் (strcmp(gotoblas_corename(), "armv9sme") == 0
#எனில் defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  எனில் (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        திரும்பு;
  }
இல்லையெனில் எனில் (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        திரும்பு;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#எனில் !defined(BFLOAT16)  && !defined(HFLOAT16)
#எனில் defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
எனில் (strcmp(gotoblas_corename(), "armv9sme") == 0
#எனில் defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
எழுத்து* TA,*TB;
      எனில் (transa & 1) 
        TA = "T";
	இல்லையெனில் 
	TA= "N";
      எனில் (transb & 1) 
        TB = "T";
	இல்லையெனில் 
	TB= "N";
#ifndef COMPLEX
      எனில் (transa == 3) 
	TA= "T";
      எனில் (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#இல்லையெனில்
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#இல்லையெனில்
      எனில் (transa == 2) 
	TA= "R";
      எனில் (transb == 2) 
	TB= "R";
      எனில் (transa == 3) 
	TA= "C";
      எனில் (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
தசம _Complex c_al={al[0],al[1]};
தசம _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#இல்லையெனில்
இரட்டை _Complex c_al={al[0],al[1]};
இரட்டை _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	திரும்பு;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#எனில் defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#எனில் defined(DYNAMIC_ARCH)
  எனில் (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    திரும்பு;
  }
#endif
#எனில் !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  எனில் (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    திரும்பு;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  எனில் ((args.m == 0) || (args.n == 0)) திரும்பு;

#எனில் 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#எனில் (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#இல்லையெனில்
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#எனில் (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#இல்லையெனில்
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#எனில் defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#எனில் defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#இல்லையெனில்
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  நிலையான மாறிலி எழுத்து gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  எனில் (args.k != 0) {
    எனில் (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      முழு trans_idx = transa | ((transb & 2) << 1);
      எழுத்து NT = gemv_trans[trans_idx];
      எனில் (transa & 1) {
        m = args.k;
        n = args.m;
      }
      எனில் (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      எனில் (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        திரும்பு;
      }
    }
    எனில் (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      முழு trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      எழுத்து NT = gemv_trans[trans_idx];
      எனில் (transa & 1) {
        inc_x = 1;
      }
      எனில் (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      எனில் (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        திரும்பு;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#எனில் USE_SMALL_MATRIX_OPT
எனில் (*(FLOAT*)args.alpha != 0.0) {
#எனில் !defined(COMPLEX)
  எனில்(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  எனில்(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }இல்லையெனில்{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  திரும்பு;
  }
#இல்லையெனில்
  எனில்(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  எனில்(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }இல்லையெனில்{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  திரும்பு;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#எனில் defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#இல்லையெனில்
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#எனில் defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (இரட்டை) args.m * (இரட்டை) args.n * (இரட்டை) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = வெற்று;

 எனில் (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, வெற்று, வெற்று, sa, sb, 0);

#ifdef SMP

  } இல்லையெனில் {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_இலக்கம்_nodes();

      எனில் ((nodes > 1) && get_முனை_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, வெற்று, வெற்று, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } இல்லையெனில் {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, வெற்று, வெற்று, sa, sb, 0);

#இல்லையெனில்

	GEMM_THREAD(mode, &args, வெற்று, வெற்று, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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

  திரும்பு;
}
