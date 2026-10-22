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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <stdlib.h>
#ಅನ್ತರ್ಭಾವ <stdbool.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#ನೋಚೇದ್ಯದಿ defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#ಅನ್ಯಥಾ
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#ನೋಚೇದ್ಯದಿ defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#ಅನ್ಯಥಾ
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#ಯದಿ !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#ಅನ್ಯಥಾ
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#ಯದಿ defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#ಅನ್ಯಥಾ
#define USE_SMALL_MATRIX_OPT 0
#endif

#ಯದಿ USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((ಶೂನ್ಯ *)(table[idx]))
#ಅನ್ಯಥಾ
#define SMALL_KERNEL_ADDR(table, idx) ((ಶೂನ್ಯ *)(*(uintptr_t *)((ವರ್ಣ *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
ಸ್ಥಿರ size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ಸ್ಥಿರ size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ಅನ್ಯಥಾ

ಸ್ಥಿರ size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ಸ್ಥಿರ size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#ಯದಿ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ openblas_amxtile_permission = 0;
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ init_amxtile_permission() {
  long status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  ಯದಿ (status != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    ಪ್ರತ್ಯಾಗಮ -1;
  }
  openblas_amxtile_permission = 1;
  ಪ್ರತ್ಯಾಗಮ 0;
}
#endif

#ifdef DYNAMIC_ARCH
ಬಾಹ್ಯ ವರ್ಣ* gotoblas_corename(ಶೂನ್ಯ);
#endif

#ifdef SMP
#ಯದಿ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ಸ್ಥಿರ ಅನ್ತಃಸ್ಥ ಪೂರ್ಣಾಙ್ಕ get_gemm_optimal_nthreads_neoversev1(ದ್ವಿಗುಣ MNK, ಪೂರ್ಣಾಙ್ಕ ncpu) {
  ಪ್ರತ್ಯಾಗಮ
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

#ಯದಿ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ಸ್ಥಿರ ಅನ್ತಃಸ್ಥ ಪೂರ್ಣಾಙ್ಕ get_gemm_optimal_nthreads_neoversev2(ದ್ವಿಗುಣ MNK, ಪೂರ್ಣಾಙ್ಕ ncpu) {
  ಪ್ರತ್ಯಾಗಮ
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

ಸ್ಥಿರ ಅನ್ತಃಸ್ಥ ಪೂರ್ಣಾಙ್ಕ get_gemm_optimal_nthreads(ದ್ವಿಗುಣ MNK) {
  ಪೂರ್ಣಾಙ್ಕ ncpu = ಅಂಕೆ_cpu_avail(3);
#ಯದಿ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ಪ್ರತ್ಯಾಗಮ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#ನೋಚೇದ್ಯದಿ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ಪ್ರತ್ಯಾಗಮ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#ನೋಚೇದ್ಯದಿ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ಯದಿ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ಪ್ರತ್ಯಾಗಮ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  ಯದಿ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ಪ್ರತ್ಯಾಗಮ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  ಯದಿ ( MNK <= (SMP_THRESHOLD_MIN  * (ದ್ವಿಗುಣ) GEMM_MULTITHREAD_THRESHOLD) ) {
    ಪ್ರತ್ಯಾಗಮ 1;
  }
  ಅನ್ಯಥಾ {
    ಯದಿ (MNK/ncpu < SMP_THRESHOLD_MIN*(ದ್ವಿಗುಣ)GEMM_MULTITHREAD_THRESHOLD) {
      ಪ್ರತ್ಯಾಗಮ MNK/(SMP_THRESHOLD_MIN*(ದ್ವಿಗುಣ)GEMM_MULTITHREAD_THRESHOLD);
    }
    ಅನ್ಯಥಾ {
      ಪ್ರತ್ಯಾಗಮ ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *TRANSA, ವರ್ಣ *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  ಪೂರ್ಣಾಙ್ಕ transa, transb, nrowa, nrowb;
  blasint info;

  ವರ್ಣ transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  ದ್ವಿಗುಣ MNK;
#ಯದಿ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#ಯದಿ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ಪೂರ್ಣಾಙ್ಕ nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (ಶೂನ್ಯ *)a;
  args.b = (ಶೂನ್ಯ *)b;
  args.c = (ಶೂನ್ಯ *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  ಯದಿ (transA == 'N') transa = 0;
  ಯದಿ (transA == 'T') transa = 1;
#ifndef COMPLEX
  ಯದಿ (transA == 'R') transa = 0;
  ಯದಿ (transA == 'C') transa = 1;
#ಅನ್ಯಥಾ
  ಯದಿ (transA == 'R') transa = 2;
  ಯದಿ (transA == 'C') transa = 3;
#endif

  ಯದಿ (transB == 'N') transb = 0;
  ಯದಿ (transB == 'T') transb = 1;
#ifndef COMPLEX
  ಯದಿ (transB == 'R') transb = 0;
  ಯದಿ (transB == 'C') transb = 1;
#ಅನ್ಯಥಾ
  ಯದಿ (transB == 'R') transb = 2;
  ಯದಿ (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  ಯದಿ (transa & 1) nrowa = args.k;
  nrowb = args.k;
  ಯದಿ (transb & 1) nrowb = args.n;

  info = 0;

  ಯದಿ (args.ldc < args.m) info = 13;
  ಯದಿ (args.ldb < nrowb)  info = 10;
  ಯದಿ (args.lda < nrowa)  info =  8;
  ಯದಿ (args.k < 0)        info =  5;
  ಯದಿ (args.n < 0)        info =  4;
  ಯದಿ (args.m < 0)        info =  3;
  ಯದಿ (transb < 0)        info =  2;
  ಯದಿ (transa < 0)        info =  1;

  ಯದಿ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_TRANSPOSE TransA, ಗಣನಾ CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint lda,
	   ಶೂನ್ಯ *vb, blasint ldb,
	   ಶೂನ್ಯ *vbeta,
	   ಶೂನ್ಯ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ಪೂರ್ಣಾಙ್ಕ transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  ದ್ವಿಗುಣ MNK;
#ಯದಿ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#ಯದಿ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ಪೂರ್ಣಾಙ್ಕ nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (ಶೂನ್ಯ *)&alpha;
  args.beta  = (ಶೂನ್ಯ *)&beta;
#ಅನ್ಯಥಾ
  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  ಯದಿ (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (ಶೂನ್ಯ *)a;
    args.b = (ಶೂನ್ಯ *)b;
    args.c = (ಶೂನ್ಯ *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    ಯದಿ (TransA == CblasNoTrans)     transa = 0;
    ಯದಿ (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    ಯದಿ (TransA == CblasConjNoTrans) transa = 0;
    ಯದಿ (TransA == CblasConjTrans)   transa = 1;
#ಅನ್ಯಥಾ
    ಯದಿ (TransA == CblasConjNoTrans) transa = 2;
    ಯದಿ (TransA == CblasConjTrans)   transa = 3;
#endif
    ಯದಿ (TransB == CblasNoTrans)     transb = 0;
    ಯದಿ (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    ಯದಿ (TransB == CblasConjNoTrans) transb = 0;
    ಯದಿ (TransB == CblasConjTrans)   transb = 1;
#ಅನ್ಯಥಾ
    ಯದಿ (TransB == CblasConjNoTrans) transb = 2;
    ಯದಿ (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    ಯದಿ (transa & 1) nrowa = args.k;
    nrowb = args.k;
    ಯದಿ (transb & 1) nrowb = args.n;

    info = -1;

    ಯದಿ (args.ldc < args.m) info = 13;
    ಯದಿ (args.ldb < nrowb)  info = 10;
    ಯದಿ (args.lda < nrowa)  info =  8;
    ಯದಿ (args.k < 0)        info =  5;
    ಯದಿ (args.n < 0)        info =  4;
    ಯದಿ (args.m < 0)        info =  3;
    ಯದಿ (transb < 0)        info =  2;
    ಯದಿ (transa < 0)        info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (ಶೂನ್ಯ *)b;
    args.b = (ಶೂನ್ಯ *)a;
    args.c = (ಶೂನ್ಯ *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    ಯದಿ (TransB == CblasNoTrans)     transa = 0;
    ಯದಿ (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    ಯದಿ (TransB == CblasConjNoTrans) transa = 0;
    ಯದಿ (TransB == CblasConjTrans)   transa = 1;
#ಅನ್ಯಥಾ
    ಯದಿ (TransB == CblasConjNoTrans) transa = 2;
    ಯದಿ (TransB == CblasConjTrans)   transa = 3;
#endif
    ಯದಿ (TransA == CblasNoTrans)     transb = 0;
    ಯದಿ (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    ಯದಿ (TransA == CblasConjNoTrans) transb = 0;
    ಯದಿ (TransA == CblasConjTrans)   transb = 1;
#ಅನ್ಯಥಾ
    ಯದಿ (TransA == CblasConjNoTrans) transb = 2;
    ಯದಿ (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    ಯದಿ (transa & 1) nrowa = args.k;
    nrowb = args.k;
    ಯದಿ (transb & 1) nrowb = args.n;

    info = -1;

    ಯದಿ (args.ldc < args.m) info = 13;
    ಯದಿ (args.ldb < nrowb)  info = 10;
    ಯದಿ (args.lda < nrowa)  info =  8;
    ಯದಿ (args.k < 0)        info =  5;
    ಯದಿ (args.n < 0)        info =  4;
    ಯದಿ (args.m < 0)        info =  3;
    ಯದಿ (transb < 0)        info =  2;
    ಯದಿ (transa < 0)        info =  1;

  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }


  ಯದಿ ((args.m == 0) || (args.n == 0)) ಪ್ರತ್ಯಾಗಮ;
#ಯದಿ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಯದಿ defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ಯದಿ defined(DYNAMIC_ARCH)
  ಯದಿ (support_avx512() )
#endif
  ಯದಿ (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
  }
#endif
#ಯದಿ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ಯದಿ defined(DYNAMIC_ARCH)
ಯದಿ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಯದಿ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  ಯದಿ (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
  }
ಅನ್ಯಥಾ ಯದಿ (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಪ್ರತ್ಯಾಗಮ;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#ಯದಿ !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಯದಿ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ಯದಿ defined(DYNAMIC_ARCH)
ಯದಿ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಯದಿ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
ವರ್ಣ* TA,*TB;
      ಯದಿ (transa & 1) 
        TA = "T";
	ಅನ್ಯಥಾ 
	TA= "N";
      ಯದಿ (transb & 1) 
        TB = "T";
	ಅನ್ಯಥಾ 
	TB= "N";
#ifndef COMPLEX
      ಯದಿ (transa == 3) 
	TA= "T";
      ಯದಿ (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#ಅನ್ಯಥಾ
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#ಅನ್ಯಥಾ
      ಯದಿ (transa == 2) 
	TA= "R";
      ಯದಿ (transb == 2) 
	TB= "R";
      ಯದಿ (transa == 3) 
	TA= "C";
      ಯದಿ (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
ಭಿನ್ನ _Complex c_al={al[0],al[1]};
ಭಿನ್ನ _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#ಅನ್ಯಥಾ
ದ್ವಿಗುಣ _Complex c_al={al[0],al[1]};
ದ್ವಿಗುಣ _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	ಪ್ರತ್ಯಾಗಮ;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#ಯದಿ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ಯದಿ defined(DYNAMIC_ARCH)
  ಯದಿ (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ಪ್ರತ್ಯಾಗಮ;
  }
#endif
#ಯದಿ !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  ಯದಿ (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ಪ್ರತ್ಯಾಗಮ;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  ಯದಿ ((args.m == 0) || (args.n == 0)) ಪ್ರತ್ಯಾಗಮ;

#ಯದಿ 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#ಯದಿ (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#ಅನ್ಯಥಾ
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#ಯದಿ (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#ಅನ್ಯಥಾ
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#ಯದಿ defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#ಯದಿ defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#ಅನ್ಯಥಾ
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ಸ್ಥಿರ ಧ್ರುವ ವರ್ಣ gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  ಯದಿ (args.k != 0) {
    ಯದಿ (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      ಪೂರ್ಣಾಙ್ಕ trans_idx = transa | ((transb & 2) << 1);
      ವರ್ಣ NT = gemv_trans[trans_idx];
      ಯದಿ (transa & 1) {
        m = args.k;
        n = args.m;
      }
      ಯದಿ (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      ಯದಿ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        ಪ್ರತ್ಯಾಗಮ;
      }
    }
    ಯದಿ (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      ಪೂರ್ಣಾಙ್ಕ trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      ವರ್ಣ NT = gemv_trans[trans_idx];
      ಯದಿ (transa & 1) {
        inc_x = 1;
      }
      ಯದಿ (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      ಯದಿ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        ಪ್ರತ್ಯಾಗಮ;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ಯದಿ USE_SMALL_MATRIX_OPT
ಯದಿ (*(FLOAT*)args.alpha != 0.0) {
#ಯದಿ !defined(COMPLEX)
  ಯದಿ(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  ಯದಿ(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }ಅನ್ಯಥಾ{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  ಪ್ರತ್ಯಾಗಮ;
  }
#ಅನ್ಯಥಾ
  ಯದಿ(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  ಯದಿ(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }ಅನ್ಯಥಾ{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  ಪ್ರತ್ಯಾಗಮ;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#ಯದಿ defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#ಅನ್ಯಥಾ
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ಯದಿ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (ದ್ವಿಗುಣ) args.m * (ದ್ವಿಗುಣ) args.n * (ದ್ವಿಗುಣ) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = ರಿಕ್ತ;

 ಯದಿ (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ifdef SMP

  } ಅನ್ಯಥಾ {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_ಅಂಕೆ_nodes();

      ಯದಿ ((nodes > 1) && get_ಗಂಟು_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, ರಿಕ್ತ, ರಿಕ್ತ, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } ಅನ್ಯಥಾ {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ಅನ್ಯಥಾ

	GEMM_THREAD(mode, &args, ರಿಕ್ತ, ರಿಕ್ತ, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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

  ಪ್ರತ್ಯಾಗಮ;
}
