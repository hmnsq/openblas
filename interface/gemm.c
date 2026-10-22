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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ <stdlib.h>
#অন্তর্ভাৱ <stdbool.h>
#অন্তর্ভাৱ "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#নোচেদ্যদি defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#অন্যথা
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#নোচেদ্যদি defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#অন্যথা
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#অন্যথা
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#অন্যথা
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#অন্যথা
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

স্থির পূর্ণাঙ্ক (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#যদি !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#অন্যথা
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#যদি defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#অন্যথা
#define USE_SMALL_MATRIX_OPT 0
#endif

#যদি USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((শূন্য *)(table[idx]))
#অন্যথা
#define SMALL_KERNEL_ADDR(table, idx) ((শূন্য *)(*(uintptr_t *)((ৱর্ণ *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
স্থির size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


স্থির size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#অন্যথা

স্থির size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

স্থির size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#যদি defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
স্থির পূর্ণাঙ্ক openblas_amxtile_permission = 0;
স্থির পূর্ণাঙ্ক init_amxtile_permission() {
  long status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  যদি (status != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    ঘূৰাই -1;
  }
  openblas_amxtile_permission = 1;
  ঘূৰাই 0;
}
#endif

#ifdef DYNAMIC_ARCH
বাহ্য ৱর্ণ* gotoblas_corename(শূন্য);
#endif

#ifdef SMP
#যদি defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
স্থির অন্তঃস্থ পূর্ণাঙ্ক get_gemm_optimal_nthreads_neoversev1(দ্ৱিগুণ MNK, পূর্ণাঙ্ক ncpu) {
  ঘূৰাই
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

#যদি defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
স্থির অন্তঃস্থ পূর্ণাঙ্ক get_gemm_optimal_nthreads_neoversev2(দ্ৱিগুণ MNK, পূর্ণাঙ্ক ncpu) {
  ঘূৰাই
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

স্থির অন্তঃস্থ পূর্ণাঙ্ক get_gemm_optimal_nthreads(দ্ৱিগুণ MNK) {
  পূর্ণাঙ্ক ncpu = অঙ্ক_cpu_avail(3);
#যদি defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ঘূৰাই get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#নোচেদ্যদি defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ঘূৰাই get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#নোচেদ্যদি defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  যদি (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ঘূৰাই get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  যদি (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ঘূৰাই get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  যদি ( MNK <= (SMP_THRESHOLD_MIN  * (দ্ৱিগুণ) GEMM_MULTITHREAD_THRESHOLD) ) {
    ঘূৰাই 1;
  }
  অন্যথা {
    যদি (MNK/ncpu < SMP_THRESHOLD_MIN*(দ্ৱিগুণ)GEMM_MULTITHREAD_THRESHOLD) {
      ঘূৰাই MNK/(SMP_THRESHOLD_MIN*(দ্ৱিগুণ)GEMM_MULTITHREAD_THRESHOLD);
    }
    অন্যথা {
      ঘূৰাই ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(ৱর্ণ *TRANSA, ৱর্ণ *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  পূর্ণাঙ্ক transa, transb, nrowa, nrowb;
  blasint info;

  ৱর্ণ transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  দ্ৱিগুণ MNK;
#যদি defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#যদি defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  পূর্ণাঙ্ক nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (শূন্য *)a;
  args.b = (শূন্য *)b;
  args.c = (শূন্য *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (শূন্য *)alpha;
  args.beta  = (শূন্য *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  যদি (transA == 'N') transa = 0;
  যদি (transA == 'T') transa = 1;
#ifndef COMPLEX
  যদি (transA == 'R') transa = 0;
  যদি (transA == 'C') transa = 1;
#অন্যথা
  যদি (transA == 'R') transa = 2;
  যদি (transA == 'C') transa = 3;
#endif

  যদি (transB == 'N') transb = 0;
  যদি (transB == 'T') transb = 1;
#ifndef COMPLEX
  যদি (transB == 'R') transb = 0;
  যদি (transB == 'C') transb = 1;
#অন্যথা
  যদি (transB == 'R') transb = 2;
  যদি (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  যদি (transa & 1) nrowa = args.k;
  nrowb = args.k;
  যদি (transb & 1) nrowb = args.n;

  info = 0;

  যদি (args.ldc < args.m) info = 13;
  যদি (args.ldb < nrowb)  info = 10;
  যদি (args.lda < nrowa)  info =  8;
  যদি (args.k < 0)        info =  5;
  যদি (args.n < 0)        info =  4;
  যদি (args.m < 0)        info =  3;
  যদি (transb < 0)        info =  2;
  যদি (transa < 0)        info =  1;

  যদি (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_TRANSPOSE TransA, গণনা CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#অন্যথা
	   শূন্য *valpha,
	   শূন্য *va, blasint lda,
	   শূন্য *vb, blasint ldb,
	   শূন্য *vbeta,
	   শূন্য *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  পূর্ণাঙ্ক transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  দ্ৱিগুণ MNK;
#যদি defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#যদি defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  পূর্ণাঙ্ক nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (শূন্য *)&alpha;
  args.beta  = (শূন্য *)&beta;
#অন্যথা
  args.alpha = (শূন্য *)alpha;
  args.beta  = (শূন্য *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  যদি (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (শূন্য *)a;
    args.b = (শূন্য *)b;
    args.c = (শূন্য *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    যদি (TransA == CblasNoTrans)     transa = 0;
    যদি (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    যদি (TransA == CblasConjNoTrans) transa = 0;
    যদি (TransA == CblasConjTrans)   transa = 1;
#অন্যথা
    যদি (TransA == CblasConjNoTrans) transa = 2;
    যদি (TransA == CblasConjTrans)   transa = 3;
#endif
    যদি (TransB == CblasNoTrans)     transb = 0;
    যদি (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    যদি (TransB == CblasConjNoTrans) transb = 0;
    যদি (TransB == CblasConjTrans)   transb = 1;
#অন্যথা
    যদি (TransB == CblasConjNoTrans) transb = 2;
    যদি (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    যদি (transa & 1) nrowa = args.k;
    nrowb = args.k;
    যদি (transb & 1) nrowb = args.n;

    info = -1;

    যদি (args.ldc < args.m) info = 13;
    যদি (args.ldb < nrowb)  info = 10;
    যদি (args.lda < nrowa)  info =  8;
    যদি (args.k < 0)        info =  5;
    যদি (args.n < 0)        info =  4;
    যদি (args.m < 0)        info =  3;
    যদি (transb < 0)        info =  2;
    যদি (transa < 0)        info =  1;
  }

  যদি (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (শূন্য *)b;
    args.b = (শূন্য *)a;
    args.c = (শূন্য *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    যদি (TransB == CblasNoTrans)     transa = 0;
    যদি (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    যদি (TransB == CblasConjNoTrans) transa = 0;
    যদি (TransB == CblasConjTrans)   transa = 1;
#অন্যথা
    যদি (TransB == CblasConjNoTrans) transa = 2;
    যদি (TransB == CblasConjTrans)   transa = 3;
#endif
    যদি (TransA == CblasNoTrans)     transb = 0;
    যদি (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    যদি (TransA == CblasConjNoTrans) transb = 0;
    যদি (TransA == CblasConjTrans)   transb = 1;
#অন্যথা
    যদি (TransA == CblasConjNoTrans) transb = 2;
    যদি (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    যদি (transa & 1) nrowa = args.k;
    nrowb = args.k;
    যদি (transb & 1) nrowb = args.n;

    info = -1;

    যদি (args.ldc < args.m) info = 13;
    যদি (args.ldb < nrowb)  info = 10;
    যদি (args.lda < nrowa)  info =  8;
    যদি (args.k < 0)        info =  5;
    যদি (args.n < 0)        info =  4;
    যদি (args.m < 0)        info =  3;
    যদি (transb < 0)        info =  2;
    যদি (transa < 0)        info =  1;

  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }


  যদি ((args.m == 0) || (args.n == 0)) ঘূৰাই;
#যদি !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#যদি defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#যদি defined(DYNAMIC_ARCH)
  যদি (support_avx512() )
#endif
  যদি (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        ঘূৰাই;
  }
#endif
#যদি defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#যদি defined(DYNAMIC_ARCH)
যদি (strcmp(gotoblas_corename(), "armv9sme") == 0
#যদি defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  যদি (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        ঘূৰাই;
  }
অন্যথা যদি (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ঘূৰাই;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#যদি !defined(BFLOAT16)  && !defined(HFLOAT16)
#যদি defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#যদি defined(DYNAMIC_ARCH)
যদি (strcmp(gotoblas_corename(), "armv9sme") == 0
#যদি defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
ৱর্ণ* TA,*TB;
      যদি (transa & 1) 
        TA = "T";
	অন্যথা 
	TA= "N";
      যদি (transb & 1) 
        TB = "T";
	অন্যথা 
	TB= "N";
#ifndef COMPLEX
      যদি (transa == 3) 
	TA= "T";
      যদি (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#অন্যথা
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#অন্যথা
      যদি (transa == 2) 
	TA= "R";
      যদি (transb == 2) 
	TB= "R";
      যদি (transa == 3) 
	TA= "C";
      যদি (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
ভিন্ন _Complex c_al={al[0],al[1]};
ভিন্ন _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#অন্যথা
দ্ৱিগুণ _Complex c_al={al[0],al[1]};
দ্ৱিগুণ _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	ঘূৰাই;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#যদি defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#যদি defined(DYNAMIC_ARCH)
  যদি (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ঘূৰাই;
  }
#endif
#যদি !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  যদি (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ঘূৰাই;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  যদি ((args.m == 0) || (args.n == 0)) ঘূৰাই;

#যদি 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#যদি (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#অন্যথা
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#যদি (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#অন্যথা
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#যদি defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#যদি defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#অন্যথা
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  স্থির ধ্রুৱ ৱর্ণ gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  যদি (args.k != 0) {
    যদি (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      পূর্ণাঙ্ক trans_idx = transa | ((transb & 2) << 1);
      ৱর্ণ NT = gemv_trans[trans_idx];
      যদি (transa & 1) {
        m = args.k;
        n = args.m;
      }
      যদি (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      যদি (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        ঘূৰাই;
      }
    }
    যদি (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      পূর্ণাঙ্ক trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      ৱর্ণ NT = gemv_trans[trans_idx];
      যদি (transa & 1) {
        inc_x = 1;
      }
      যদি (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      যদি (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        ঘূৰাই;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#যদি USE_SMALL_MATRIX_OPT
যদি (*(FLOAT*)args.alpha != 0.0) {
#যদি !defined(COMPLEX)
  যদি(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  যদি(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }অন্যথা{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  ঘূৰাই;
  }
#অন্যথা
  যদি(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  যদি(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }অন্যথা{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  ঘূৰাই;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#যদি defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#অন্যথা
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#যদি defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (দ্ৱিগুণ) args.m * (দ্ৱিগুণ) args.n * (দ্ৱিগুণ) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = রিক্ত;

 যদি (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, রিক্ত, রিক্ত, sa, sb, 0);

#ifdef SMP

  } অন্যথা {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_অঙ্ক_nodes();

      যদি ((nodes > 1) && get_গ্ৰন্থি_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, রিক্ত, রিক্ত, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } অন্যথা {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, রিক্ত, রিক্ত, sa, sb, 0);

#অন্যথা

	GEMM_THREAD(mode, &args, রিক্ত, রিক্ত, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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

  ঘূৰাই;
}
