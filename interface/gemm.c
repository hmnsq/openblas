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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ <stdbool.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#ନଚେତଯଦି defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#ଅନ୍ୟଥା
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#ନଚେତଯଦି defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#ଅନ୍ୟଥା
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#ଅନ୍ୟଥା
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#ଯଦି !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#ଯଦି defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#ଅନ୍ୟଥା
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#ଯଦି defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#ଯଦି defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#ଅନ୍ୟଥା
#define USE_SMALL_MATRIX_OPT 0
#endif

#ଯଦି USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((ଶୂନ୍ୟ *)(table[idx]))
#ଅନ୍ୟଥା
#define SMALL_KERNEL_ADDR(table, idx) ((ଶୂନ୍ୟ *)(*(uintptr_t *)((ଅକ୍ଷର *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
ସ୍ଥିର size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ସ୍ଥିର size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ଅନ୍ୟଥା

ସ୍ଥିର size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ସ୍ଥିର size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#ଯଦି defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା openblas_amxtile_permission = 0;
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା init_amxtile_permission() {
  long status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  ଯଦି (status != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    ଫେରନ୍ତୁ -1;
  }
  openblas_amxtile_permission = 1;
  ଫେରନ୍ତୁ 0;
}
#endif

#ifdef DYNAMIC_ARCH
ବାହ୍ୟ ଅକ୍ଷର* gotoblas_corename(ଶୂନ୍ୟ);
#endif

#ifdef SMP
#ଯଦି defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ସ୍ଥିର ଧାଡ଼ିଭିତରେ ପୂର୍ଣ୍ଣସଂଖ୍ୟା get_gemm_optimal_nthreads_neoversev1(ଦ୍ୱିଗୁଣ MNK, ପୂର୍ଣ୍ଣସଂଖ୍ୟା ncpu) {
  ଫେରନ୍ତୁ
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

#ଯଦି defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ସ୍ଥିର ଧାଡ଼ିଭିତରେ ପୂର୍ଣ୍ଣସଂଖ୍ୟା get_gemm_optimal_nthreads_neoversev2(ଦ୍ୱିଗୁଣ MNK, ପୂର୍ଣ୍ଣସଂଖ୍ୟା ncpu) {
  ଫେରନ୍ତୁ
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

ସ୍ଥିର ଧାଡ଼ିଭିତରେ ପୂର୍ଣ୍ଣସଂଖ୍ୟା get_gemm_optimal_nthreads(ଦ୍ୱିଗୁଣ MNK) {
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା ncpu = ଅଙ୍କ_cpu_avail(3);
#ଯଦି defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ଫେରନ୍ତୁ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#ନଚେତଯଦି defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ଫେରନ୍ତୁ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#ନଚେତଯଦି defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ଯଦି (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ଫେରନ୍ତୁ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  ଯଦି (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ଫେରନ୍ତୁ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  ଯଦି ( MNK <= (SMP_THRESHOLD_MIN  * (ଦ୍ୱିଗୁଣ) GEMM_MULTITHREAD_THRESHOLD) ) {
    ଫେରନ୍ତୁ 1;
  }
  ଅନ୍ୟଥା {
    ଯଦି (MNK/ncpu < SMP_THRESHOLD_MIN*(ଦ୍ୱିଗୁଣ)GEMM_MULTITHREAD_THRESHOLD) {
      ଫେରନ୍ତୁ MNK/(SMP_THRESHOLD_MIN*(ଦ୍ୱିଗୁଣ)GEMM_MULTITHREAD_THRESHOLD);
    }
    ଅନ୍ୟଥା {
      ଫେରନ୍ତୁ ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *TRANSA, ଅକ୍ଷର *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା transa, transb, nrowa, nrowb;
  blasint info;

  ଅକ୍ଷର transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  ଦ୍ୱିଗୁଣ MNK;
#ଯଦି defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#ଯଦି defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (ଶୂନ୍ୟ *)a;
  args.b = (ଶୂନ୍ୟ *)b;
  args.c = (ଶୂନ୍ୟ *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (ଶୂନ୍ୟ *)alpha;
  args.beta  = (ଶୂନ୍ୟ *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  ଯଦି (transA == 'N') transa = 0;
  ଯଦି (transA == 'T') transa = 1;
#ifndef COMPLEX
  ଯଦି (transA == 'R') transa = 0;
  ଯଦି (transA == 'C') transa = 1;
#ଅନ୍ୟଥା
  ଯଦି (transA == 'R') transa = 2;
  ଯଦି (transA == 'C') transa = 3;
#endif

  ଯଦି (transB == 'N') transb = 0;
  ଯଦି (transB == 'T') transb = 1;
#ifndef COMPLEX
  ଯଦି (transB == 'R') transb = 0;
  ଯଦି (transB == 'C') transb = 1;
#ଅନ୍ୟଥା
  ଯଦି (transB == 'R') transb = 2;
  ଯଦି (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  ଯଦି (transa & 1) nrowa = args.k;
  nrowb = args.k;
  ଯଦି (transb & 1) nrowb = args.n;

  info = 0;

  ଯଦି (args.ldc < args.m) info = 13;
  ଯଦି (args.ldb < nrowb)  info = 10;
  ଯଦି (args.lda < nrowa)  info =  8;
  ଯଦି (args.k < 0)        info =  5;
  ଯଦି (args.n < 0)        info =  4;
  ଯଦି (args.m < 0)        info =  3;
  ଯଦି (transb < 0)        info =  2;
  ଯଦି (transa < 0)        info =  1;

  ଯଦି (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_TRANSPOSE TransA, ଗଣନା CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *valpha,
	   ଶୂନ୍ୟ *va, blasint lda,
	   ଶୂନ୍ୟ *vb, blasint ldb,
	   ଶୂନ୍ୟ *vbeta,
	   ଶୂନ୍ୟ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  ଦ୍ୱିଗୁଣ MNK;
#ଯଦି defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#ଯଦି defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (ଶୂନ୍ୟ *)&alpha;
  args.beta  = (ଶୂନ୍ୟ *)&beta;
#ଅନ୍ୟଥା
  args.alpha = (ଶୂନ୍ୟ *)alpha;
  args.beta  = (ଶୂନ୍ୟ *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  ଯଦି (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (ଶୂନ୍ୟ *)a;
    args.b = (ଶୂନ୍ୟ *)b;
    args.c = (ଶୂନ୍ୟ *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    ଯଦି (TransA == CblasNoTrans)     transa = 0;
    ଯଦି (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    ଯଦି (TransA == CblasConjNoTrans) transa = 0;
    ଯଦି (TransA == CblasConjTrans)   transa = 1;
#ଅନ୍ୟଥା
    ଯଦି (TransA == CblasConjNoTrans) transa = 2;
    ଯଦି (TransA == CblasConjTrans)   transa = 3;
#endif
    ଯଦି (TransB == CblasNoTrans)     transb = 0;
    ଯଦି (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    ଯଦି (TransB == CblasConjNoTrans) transb = 0;
    ଯଦି (TransB == CblasConjTrans)   transb = 1;
#ଅନ୍ୟଥା
    ଯଦି (TransB == CblasConjNoTrans) transb = 2;
    ଯଦି (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    ଯଦି (transa & 1) nrowa = args.k;
    nrowb = args.k;
    ଯଦି (transb & 1) nrowb = args.n;

    info = -1;

    ଯଦି (args.ldc < args.m) info = 13;
    ଯଦି (args.ldb < nrowb)  info = 10;
    ଯଦି (args.lda < nrowa)  info =  8;
    ଯଦି (args.k < 0)        info =  5;
    ଯଦି (args.n < 0)        info =  4;
    ଯଦି (args.m < 0)        info =  3;
    ଯଦି (transb < 0)        info =  2;
    ଯଦି (transa < 0)        info =  1;
  }

  ଯଦି (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (ଶୂନ୍ୟ *)b;
    args.b = (ଶୂନ୍ୟ *)a;
    args.c = (ଶୂନ୍ୟ *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    ଯଦି (TransB == CblasNoTrans)     transa = 0;
    ଯଦି (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    ଯଦି (TransB == CblasConjNoTrans) transa = 0;
    ଯଦି (TransB == CblasConjTrans)   transa = 1;
#ଅନ୍ୟଥା
    ଯଦି (TransB == CblasConjNoTrans) transa = 2;
    ଯଦି (TransB == CblasConjTrans)   transa = 3;
#endif
    ଯଦି (TransA == CblasNoTrans)     transb = 0;
    ଯଦି (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    ଯଦି (TransA == CblasConjNoTrans) transb = 0;
    ଯଦି (TransA == CblasConjTrans)   transb = 1;
#ଅନ୍ୟଥା
    ଯଦି (TransA == CblasConjNoTrans) transb = 2;
    ଯଦି (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    ଯଦି (transa & 1) nrowa = args.k;
    nrowb = args.k;
    ଯଦି (transb & 1) nrowb = args.n;

    info = -1;

    ଯଦି (args.ldc < args.m) info = 13;
    ଯଦି (args.ldb < nrowb)  info = 10;
    ଯଦି (args.lda < nrowa)  info =  8;
    ଯଦି (args.k < 0)        info =  5;
    ଯଦି (args.n < 0)        info =  4;
    ଯଦି (args.m < 0)        info =  3;
    ଯଦି (transb < 0)        info =  2;
    ଯଦି (transa < 0)        info =  1;

  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }


  ଯଦି ((args.m == 0) || (args.n == 0)) ଫେରନ୍ତୁ;
#ଯଦି !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ଯଦି defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ଯଦି defined(DYNAMIC_ARCH)
  ଯଦି (support_avx512() )
#endif
  ଯଦି (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        ଫେରନ୍ତୁ;
  }
#endif
#ଯଦି defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ଯଦି defined(DYNAMIC_ARCH)
ଯଦି (strcmp(gotoblas_corename(), "armv9sme") == 0
#ଯଦି defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  ଯଦି (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        ଫେରନ୍ତୁ;
  }
ଅନ୍ୟଥା ଯଦି (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ଫେରନ୍ତୁ;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#ଯଦି !defined(BFLOAT16)  && !defined(HFLOAT16)
#ଯଦି defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ଯଦି defined(DYNAMIC_ARCH)
ଯଦି (strcmp(gotoblas_corename(), "armv9sme") == 0
#ଯଦି defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
ଅକ୍ଷର* TA,*TB;
      ଯଦି (transa & 1) 
        TA = "T";
	ଅନ୍ୟଥା 
	TA= "N";
      ଯଦି (transb & 1) 
        TB = "T";
	ଅନ୍ୟଥା 
	TB= "N";
#ifndef COMPLEX
      ଯଦି (transa == 3) 
	TA= "T";
      ଯଦି (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#ଅନ୍ୟଥା
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#ଅନ୍ୟଥା
      ଯଦି (transa == 2) 
	TA= "R";
      ଯଦି (transb == 2) 
	TB= "R";
      ଯଦି (transa == 3) 
	TA= "C";
      ଯଦି (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
ଦଶମିକ _Complex c_al={al[0],al[1]};
ଦଶମିକ _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#ଅନ୍ୟଥା
ଦ୍ୱିଗୁଣ _Complex c_al={al[0],al[1]};
ଦ୍ୱିଗୁଣ _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	ଫେରନ୍ତୁ;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#ଯଦି defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ଯଦି defined(DYNAMIC_ARCH)
  ଯଦି (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ଫେରନ୍ତୁ;
  }
#endif
#ଯଦି !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  ଯଦି (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ଫେରନ୍ତୁ;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  ଯଦି ((args.m == 0) || (args.n == 0)) ଫେରନ୍ତୁ;

#ଯଦି 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#ଯଦି (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#ଅନ୍ୟଥା
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#ଯଦି (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#ଅନ୍ୟଥା
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#ଯଦି defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#ଯଦି defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#ଅନ୍ୟଥା
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ସ୍ଥିର ଅଚଳ ଅକ୍ଷର gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  ଯଦି (args.k != 0) {
    ଯଦି (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans_idx = transa | ((transb & 2) << 1);
      ଅକ୍ଷର NT = gemv_trans[trans_idx];
      ଯଦି (transa & 1) {
        m = args.k;
        n = args.m;
      }
      ଯଦି (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      ଯଦି (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        ଫେରନ୍ତୁ;
      }
    }
    ଯଦି (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      ଅକ୍ଷର NT = gemv_trans[trans_idx];
      ଯଦି (transa & 1) {
        inc_x = 1;
      }
      ଯଦି (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      ଯଦି (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        ଫେରନ୍ତୁ;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ଯଦି USE_SMALL_MATRIX_OPT
ଯଦି (*(FLOAT*)args.alpha != 0.0) {
#ଯଦି !defined(COMPLEX)
  ଯଦି(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  ଯଦି(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }ଅନ୍ୟଥା{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  ଫେରନ୍ତୁ;
  }
#ଅନ୍ୟଥା
  ଯଦି(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  ଯଦି(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }ଅନ୍ୟଥା{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  ଫେରନ୍ତୁ;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#ଯଦି defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#ଅନ୍ୟଥା
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ଯଦି defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (ଦ୍ୱିଗୁଣ) args.m * (ଦ୍ୱିଗୁଣ) args.n * (ଦ୍ୱିଗୁଣ) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = ଖାଲି;

 ଯଦି (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP

  } ଅନ୍ୟଥା {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_ଅଙ୍କ_nodes();

      ଯଦି ((nodes > 1) && get_ଗଣ୍ଠି_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, ଖାଲି, ଖାଲି, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } ଅନ୍ୟଥା {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ଅନ୍ୟଥା

	GEMM_THREAD(mode, &args, ଖାଲି, ଖାଲି, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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

  ଫେରନ୍ତୁ;
}
