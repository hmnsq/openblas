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

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل <stdbool.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_غیر_تعریف_شدہ COMPLEX
#تعریف_کرو SMP_حد_MIN 65536.0
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMM "
#تعریف_کرو GEMV BLASFUNC(qgemv)
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMM "
#تعریف_کرو GEMV BLASFUNC(dgemv)
#وگرنہ_اگر defined(BFLOAT16)
#اگر_تعریف_شدہ BGEMM
#تعریف_کرو ERROR_NAME "BGEMM "
#تعریف_کرو GEMV BLASFUNC(bgemv)
#ورنہ
#تعریف_کرو ERROR_NAME "SBGEMM "
#تعریف_کرو GEMV BLASFUNC(sbgemv)
#اگر_ختم
#وگرنہ_اگر defined(HFLOAT16)
#تعریف_کرو ERROR_NAME "SHGEMM "
#ورنہ
#تعریف_کرو ERROR_NAME "SGEMM "
#تعریف_کرو GEMV BLASFUNC(sgemv)
#اگر_ختم
#ورنہ
#تعریف_کرو SMP_حد_MIN 8192.0
#اگر_غیر_تعریف_شدہ GEMM3M
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM "
#تعریف_کرو GEMV BLASFUNC(xgemv)
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM "
#تعریف_کرو GEMV BLASFUNC(zgemv)
#ورنہ
#تعریف_کرو ERROR_NAME "CGEMM "
#تعریف_کرو GEMV BLASFUNC(cgemv)
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM3M "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM3M "
#ورنہ
#تعریف_کرو ERROR_NAME "CGEMM3M "
#اگر_ختم
#اگر_ختم
#اگر_ختم

#اگر_غیر_تعریف_شدہ GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#اگر_ختم

ساکن صحیح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#اگر !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#اگر_ختم
#ورنہ
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#اگر_ختم
#اگر_ختم
};

#اگر defined(SMALL_میٹرکس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#تعریف_کرو USE_SMALL_میٹرکس_OPT 1
#ورنہ
#تعریف_کرو USE_SMALL_میٹرکس_OPT 0
#اگر_ختم

#اگر USE_SMALL_میٹرکس_OPT
#اگر_غیر_تعریف_شدہ DYNAMIC_ARCH
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((خالی *)(جدول[idx]))
#ورنہ
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((خالی *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#اگر_ختم


#اگر_غیر_تعریف_شدہ COMPLEX
ساکن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساکن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کرو GEMM_SMALL_KERNEL_B0(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کرو GEMM_SMALL_KERNEL(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ورنہ

ساکن size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ساکن size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#تعریف_کرو ZGEMM_SMALL_KERNEL(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کرو ZGEMM_SMALL_KERNEL_B0(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#اگر_ختم
#اگر_ختم

#اگر defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#تعریف_کرو XFEATURE_XTILEDATA 18
#تعریف_کرو ARCH_REQ_XCOMP_PERM 0x1023
ساکن صحیح openblas_amxtile_permission = 0;
ساکن صحیح init_amxtile_permission() {
  طویل حالت_2 =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  اگر (حالت_2 != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    واپس -1;
  }
  openblas_amxtile_permission = 1;
  واپس 0;
}
#اگر_ختم

#اگر_تعریف_شدہ DYNAMIC_ARCH
خارجی حرف* gotoblas_corename(خالی);
#اگر_ختم

#اگر_تعریف_شدہ SMP
#اگر defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ساکن اندرونی صحیح get_gemm_optimal_nthreads_neoversev1(دوہرا MNK, صحیح ncpu) {
  واپس
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
#اگر_ختم

#اگر defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ساکن اندرونی صحیح get_gemm_optimal_nthreads_neoversev2(دوہرا MNK, صحیح ncpu) {
  واپس
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
#اگر_ختم

ساکن اندرونی صحیح get_gemm_optimal_nthreads(دوہرا MNK) {
  صحیح ncpu = عدد_cpu_avail(3);
#اگر defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  واپس get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#وگرنہ_اگر defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  واپس get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#وگرنہ_اگر defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  اگر (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    واپس get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  اگر (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    واپس get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#اگر_ختم
  اگر ( MNK <= (SMP_حد_MIN  * (دوہرا) GEMM_MULTITHREAD_حد) ) {
    واپس 1;
  }
  ورنہ {
    اگر (MNK/ncpu < SMP_حد_MIN*(دوہرا)GEMM_MULTITHREAD_حد) {
      واپس MNK/(SMP_حد_MIN*(دوہرا)GEMM_MULTITHREAD_حد);
    }
    ورنہ {
      واپس ncpu;
    }
  }
}
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *TRANSA, حرف *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t دلائل;

  صحیح transa, transb, nrowa, nrowb;
  blasint info;

  حرف transA, transB;
  IFLOAT *بفر;
  IFLOAT *sa, *sb;

#اگر_تعریف_شدہ SMP
  دوہرا MNK;
#اگر defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم
#اگر_ختم

#اگر defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  صحیح nodes;
#اگر_ختم

  PRINT_DEBUG_NAME;

  دلائل.m = *M;
  دلائل.n = *N;
  دلائل.k = *K;

  دلائل.a = (خالی *)a;
  دلائل.b = (خالی *)b;
  دلائل.c = (خالی *)c;

  دلائل.lda = *ldA;
  دلائل.ldb = *ldB;
  دلائل.ldc = *ldC;

  دلائل.alpha = (خالی *)alpha;
  دلائل.beta  = (خالی *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  اگر (transA == 'N') transa = 0;
  اگر (transA == 'T') transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (transA == 'R') transa = 0;
  اگر (transA == 'C') transa = 1;
#ورنہ
  اگر (transA == 'R') transa = 2;
  اگر (transA == 'C') transa = 3;
#اگر_ختم

  اگر (transB == 'N') transb = 0;
  اگر (transB == 'T') transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (transB == 'R') transb = 0;
  اگر (transB == 'C') transb = 1;
#ورنہ
  اگر (transB == 'R') transb = 2;
  اگر (transB == 'C') transb = 3;
#اگر_ختم
  nrowa = دلائل.m;
  اگر (transa & 1) nrowa = دلائل.k;
  nrowb = دلائل.k;
  اگر (transb & 1) nrowb = دلائل.n;

  info = 0;

  اگر (دلائل.ldc < دلائل.m) info = 13;
  اگر (دلائل.ldb < nrowb)  info = 10;
  اگر (دلائل.lda < nrowa)  info =  8;
  اگر (دلائل.k < 0)        info =  5;
  اگر (دلائل.n < 0)        info =  4;
  اگر (دلائل.m < 0)        info =  3;
  اگر (transb < 0)        info =  2;
  اگر (transa < 0)        info =  1;

  اگر (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_TRANSPOSE TransA, شمار CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#اگر_غیر_تعریف_شدہ COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ورنہ
	   خالی *valpha,
	   خالی *va, blasint lda,
	   خالی *vb, blasint ldb,
	   خالی *vbeta,
	   خالی *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#اگر_ختم

  blas_arg_t دلائل;
  صحیح transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#اگر_تعریف_شدہ SMP
  دوہرا MNK;
#اگر defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم
#اگر_ختم

#اگر defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  صحیح nodes;
#اگر_ختم

  PRINT_DEBUG_CNAME;

#اگر_غیر_تعریف_شدہ COMPLEX
  دلائل.alpha = (خالی *)&alpha;
  دلائل.beta  = (خالی *)&beta;
#ورنہ
  دلائل.alpha = (خالی *)alpha;
  دلائل.beta  = (خالی *)beta;
#اگر_ختم

  transa = -1;
  transb = -1;
  info   =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    دلائل.m = m;
    دلائل.n = n;
    دلائل.k = k;

    دلائل.a = (خالی *)a;
    دلائل.b = (خالی *)b;
    دلائل.c = (خالی *)c;

    دلائل.lda = lda;
    دلائل.ldb = ldb;
    دلائل.ldc = ldc;

    اگر (TransA == CblasNoTrans)     transa = 0;
    اگر (TransA == CblasTrans)       transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (TransA == CblasConjNoTrans) transa = 0;
    اگر (TransA == CblasConjTrans)   transa = 1;
#ورنہ
    اگر (TransA == CblasConjNoTrans) transa = 2;
    اگر (TransA == CblasConjTrans)   transa = 3;
#اگر_ختم
    اگر (TransB == CblasNoTrans)     transb = 0;
    اگر (TransB == CblasTrans)       transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (TransB == CblasConjNoTrans) transb = 0;
    اگر (TransB == CblasConjTrans)   transb = 1;
#ورنہ
    اگر (TransB == CblasConjNoTrans) transb = 2;
    اگر (TransB == CblasConjTrans)   transb = 3;
#اگر_ختم

    nrowa = دلائل.m;
    اگر (transa & 1) nrowa = دلائل.k;
    nrowb = دلائل.k;
    اگر (transb & 1) nrowb = دلائل.n;

    info = -1;

    اگر (دلائل.ldc < دلائل.m) info = 13;
    اگر (دلائل.ldb < nrowb)  info = 10;
    اگر (دلائل.lda < nrowa)  info =  8;
    اگر (دلائل.k < 0)        info =  5;
    اگر (دلائل.n < 0)        info =  4;
    اگر (دلائل.m < 0)        info =  3;
    اگر (transb < 0)        info =  2;
    اگر (transa < 0)        info =  1;
  }

  اگر (ترتیب == Cblas_صف_Major) {
    دلائل.m = n;
    دلائل.n = m;
    دلائل.k = k;

    دلائل.a = (خالی *)b;
    دلائل.b = (خالی *)a;
    دلائل.c = (خالی *)c;

    دلائل.lda = ldb;
    دلائل.ldb = lda;
    دلائل.ldc = ldc;

    اگر (TransB == CblasNoTrans)     transa = 0;
    اگر (TransB == CblasTrans)       transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (TransB == CblasConjNoTrans) transa = 0;
    اگر (TransB == CblasConjTrans)   transa = 1;
#ورنہ
    اگر (TransB == CblasConjNoTrans) transa = 2;
    اگر (TransB == CblasConjTrans)   transa = 3;
#اگر_ختم
    اگر (TransA == CblasNoTrans)     transb = 0;
    اگر (TransA == CblasTrans)       transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (TransA == CblasConjNoTrans) transb = 0;
    اگر (TransA == CblasConjTrans)   transb = 1;
#ورنہ
    اگر (TransA == CblasConjNoTrans) transb = 2;
    اگر (TransA == CblasConjTrans)   transb = 3;
#اگر_ختم

    nrowa = دلائل.m;
    اگر (transa & 1) nrowa = دلائل.k;
    nrowb = دلائل.k;
    اگر (transb & 1) nrowb = دلائل.n;

    info = -1;

    اگر (دلائل.ldc < دلائل.m) info = 13;
    اگر (دلائل.ldb < nrowb)  info = 10;
    اگر (دلائل.lda < nrowa)  info =  8;
    اگر (دلائل.k < 0)        info =  5;
    اگر (دلائل.n < 0)        info =  4;
    اگر (دلائل.m < 0)        info =  3;
    اگر (transb < 0)        info =  2;
    اگر (transa < 0)        info =  1;

  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }


  اگر ((دلائل.m == 0) || (دلائل.n == 0)) واپس;
#اگر !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#اگر defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_براہ_راست)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
  اگر (support_avx512() )
#اگر_ختم
  اگر (ترتیب == Cblas_صف_Major && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_براہ_راست_PERFORMANT(m,n,k)) {
        SGEMM_براہ_راست(m, n, k, a, lda, b, ldb, c, ldc);
        واپس;
  }
#اگر_ختم
#اگر defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_براہ_راست)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
اگر (strcmp(gotoblas_corename(), "armv9sme") == 0
#اگر defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#اگر_ختم
)
#اگر_ختم

  اگر (ترتیب == Cblas_صف_Major && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_براہ_راست_PERFORMANT(m,n,k)) {
        SGEMM_براہ_راست(m, n, k, a, lda, b, ldb, c, ldc);
        واپس;
  }
ورنہ اگر (ترتیب == Cblas_صف_Major && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_براہ_راست_PERFORMANT(m,n,k)) {
        SGEMM_براہ_راست_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        واپس;
  }

#اگر_ختم //defined arm64
#اگر_ختم //defined complex


#اگر_ختم //defined CBLAS

#اگر !defined(BFLOAT16)  && !defined(HFLOAT16)
#اگر defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_براہ_راست)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
اگر (strcmp(gotoblas_corename(), "armv9sme") == 0
#اگر defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#اگر_ختم
)
#اگر_ختم //defined dynarch
{
حرف* TA,*TB;
      اگر (transa & 1) 
        TA = "T";
	ورنہ 
	TA= "N";
      اگر (transb & 1) 
        TB = "T";
	ورنہ 
	TB= "N";
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transa == 3) 
	TA= "T";
      اگر (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)دلائل.alpha;
FLOAT* be=(FLOAT*)دلائل.beta;
#اگر_غیر_تعریف_شدہ DOUBLE
        SME_SGEMM_KERNEL(TA,TB, دلائل.m, دلائل.n, دلائل.k, al, دلائل.a, دلائل.lda, دلائل.b, دلائل.ldb, be, دلائل.c, دلائل.ldc);
#ورنہ
        SME_DGEMM_KERNEL(TA,TB, دلائل.m, دلائل.n, دلائل.k, al, دلائل.a, دلائل.lda, دلائل.b, دلائل.ldb, be, دلائل.c, دلائل.ldc);
#اگر_ختم
#ورنہ
      اگر (transa == 2) 
	TA= "R";
      اگر (transb == 2) 
	TB= "R";
      اگر (transa == 3) 
	TA= "C";
      اگر (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)دلائل.alpha;
FLOAT* be=(FLOAT*)دلائل.beta;
#اگر_غیر_تعریف_شدہ DOUBLE
اعشاری _Complex c_al={al[0],al[1]};
اعشاری _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, دلائل.m, دلائل.n, دلائل.k, al[0],al[1], دلائل.a, دلائل.lda, دلائل.b, دلائل.ldb, be[0],be[1], دلائل.c, دلائل.ldc);
#ورنہ
دوہرا _Complex c_al={al[0],al[1]};
دوہرا _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, دلائل.m, دلائل.n, دلائل.k, al[0],al[1], دلائل.a, دلائل.lda, دلائل.b, دلائل.ldb, be[0],be[1], دلائل.c, دلائل.ldc);
#اگر_ختم
#اگر_ختم
	واپس;
}
#اگر_ختم //defined arm64

#اگر_ختم //defined b/hfloat16 

#اگر defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#اگر defined(DYNAMIC_ARCH)
  اگر (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    واپس;
  }
#اگر_ختم
#اگر !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  اگر (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    واپس;
  }
#اگر_ختم
#اگر_ختم  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  اگر ((دلائل.m == 0) || (دلائل.n == 0)) واپس;

#اگر 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 دلائل.m, دلائل.n, دلائل.k, دلائل.lda, دلائل.ldb, دلائل.ldc);
#اگر_ختم

#اگر (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#تعریف_کرو BFLOAT16_GEMM_GEMV_FORWARD 1
#ورنہ
#تعریف_کرو BFLOAT16_GEMM_GEMV_FORWARD 0
#اگر_ختم
#اگر (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#تعریف_کرو HFLOAT16_GEMM_GEMV_FORWARD 1
#ورنہ
#تعریف_کرو HFLOAT16_GEMM_GEMV_FORWARD 0
#اگر_ختم

#اگر defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#اگر defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#ورنہ
  bool have_tuned_gemv = true;
#اگر_ختم
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ساکن ثابت حرف gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  اگر (دلائل.k != 0) {
    اگر (دلائل.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = دلائل.m;
      blasint n = دلائل.k;
      blasint lda = دلائل.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      صحیح trans_idx = transa | ((transb & 2) << 1);
      حرف NT = gemv_trans[trans_idx];
      اگر (transa & 1) {
        m = دلائل.k;
        n = دلائل.m;
      }
      اگر (transb & 1) {
        inc_x = دلائل.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      اگر (is_efficient_gemv) {
        GEMV(&NT, &m, &n, دلائل.alpha, دلائل.a, &lda, دلائل.b, &inc_x, دلائل.beta, دلائل.c, &inc_y);
        واپس;
      }
    }
    اگر (دلائل.m == 1) {
      blasint inc_x = دلائل.lda;
      blasint inc_y = دلائل.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = دلائل.k;
      blasint n = دلائل.n;
      blasint ldb = دلائل.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      صحیح trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      حرف NT = gemv_trans[trans_idx];
      اگر (transa & 1) {
        inc_x = 1;
      }
      اگر (transb & 1) {
        m = دلائل.n;
        n = دلائل.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      اگر (is_efficient_gemv) {
        GEMV(&NT, &m, &n, دلائل.alpha, دلائل.b, &ldb, دلائل.a, &inc_x, دلائل.beta, دلائل.c, &inc_y);
        واپس;
      }
    }
  }
#اگر_ختم

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

#اگر USE_SMALL_میٹرکس_OPT
اگر (*(FLOAT*)دلائل.alpha != 0.0) {
#اگر !defined(COMPLEX)
  اگر(GEMM_SMALL_میٹرکس_PERMIT(transa, transb, دلائل.m, دلائل.n, دلائل.k, *(FLOAT *)(دلائل.alpha), *(FLOAT *)(دلائل.beta))){

	  اگر(*(FLOAT *)(دلائل.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(دلائل.m, دلائل.n, دلائل.k, دلائل.a, دلائل.lda, *(FLOAT *)(دلائل.alpha), دلائل.b, دلائل.ldb, دلائل.c, دلائل.ldc);
	  }ورنہ{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(دلائل.m, دلائل.n, دلائل.k, دلائل.a, دلائل.lda, *(FLOAT *)(دلائل.alpha), دلائل.b, دلائل.ldb, *(FLOAT *)(دلائل.beta), دلائل.c, دلائل.ldc);
	  }
	  واپس;
  }
#ورنہ
  اگر(GEMM_SMALL_میٹرکس_PERMIT(transa, transb, دلائل.m, دلائل.n, دلائل.k, alpha[0], alpha[1], beta[0], beta[1])){
	  اگر(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(دلائل.m, دلائل.n, دلائل.k, دلائل.a, دلائل.lda, alpha[0], alpha[1], دلائل.b, دلائل.ldb, دلائل.c, دلائل.ldc);
	  }ورنہ{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(دلائل.m, دلائل.n, دلائل.k, دلائل.a, دلائل.lda, alpha[0], alpha[1], دلائل.b, دلائل.ldb, beta[0], beta[1], دلائل.c, دلائل.ldc);
	  }
	  واپس;
  }
#اگر_ختم
}
#اگر_ختم

  بفر = (XFLOAT *)blas_حافظہ_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#اگر defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)بفر + (WhereAmI() & 0xf) * GEMM_انحراف_A);
#ورنہ
  sa = (XFLOAT *)((BLASLONG)بفر +GEMM_انحراف_A);
#اگر_ختم
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);

#اگر_تعریف_شدہ SMP
#اگر defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  منوال |= (transa << BLAS_TRANSA_سرکاؤ);
  منوال |= (transb << BLAS_TRANSB_سرکاؤ);
#اگر_ختم

  MNK = (دوہرا) دلائل.m * (دوہرا) دلائل.n * (دوہرا) دلائل.k;
  دلائل.nthreads = get_gemm_optimal_nthreads(MNK);

  دلائل.common = تہی;

 اگر (دلائل.nthreads == 1) {
#اگر_ختم

    (gemm[(transb << 2) | transa])(&دلائل, تہی, تہی, sa, sb, 0);

#اگر_تعریف_شدہ SMP

  } ورنہ {

#اگر_غیر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3

#اگر_غیر_تعریف_شدہ NO_AFFINITY
      nodes = get_عدد_nodes();

      اگر ((nodes > 1) && get_گرہ_equal()) {

	دلائل.nthreads /= nodes;

	gemm_thread_mn(منوال, &دلائل, تہی, تہی, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } ورنہ {
#اگر_ختم

	(gemm[16 | (transb << 2) | transa])(&دلائل, تہی, تہی, sa, sb, 0);

#ورنہ

	GEMM_THREAD(منوال, &دلائل, تہی, تہی, gemm[(transb << 2) | transa], sa, sb, دلائل.nthreads);

#اگر_ختم

#اگر_غیر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
#اگر_غیر_تعریف_شدہ NO_AFFINITY
      }
#اگر_ختم
#اگر_ختم

#اگر_ختم

#اگر_تعریف_شدہ SMP
  }
#اگر_ختم

 blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE, دلائل.m * دلائل.k + دلائل.k * دلائل.n + دلائل.m * دلائل.n, 2 * دلائل.m * دلائل.n * دلائل.k);

  IDEBUG_انتہا;

  واپس;
}
