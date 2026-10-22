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

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن <stdbool.h>
#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_بنا_تعریف COMPLEX
#تعریف_کرو SMP_حد_MIN 65536.0
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMM "
#تعریف_کرو GEMV BLASFUNC(qgemv)
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMM "
#تعریف_کرو GEMV BLASFUNC(dgemv)
#وإلاإذا defined(BFLOAT16)
#جے_تعریف_شدہ BGEMM
#تعریف_کرو ERROR_NAME "BGEMM "
#تعریف_کرو GEMV BLASFUNC(bgemv)
#وإلا
#تعریف_کرو ERROR_NAME "SBGEMM "
#تعریف_کرو GEMV BLASFUNC(sbgemv)
#جے_انت
#وإلاإذا defined(HFLOAT16)
#تعریف_کرو ERROR_NAME "SHGEMM "
#وإلا
#تعریف_کرو ERROR_NAME "SGEMM "
#تعریف_کرو GEMV BLASFUNC(sgemv)
#جے_انت
#وإلا
#تعریف_کرو SMP_حد_MIN 8192.0
#جے_بنا_تعریف GEMM3M
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM "
#تعریف_کرو GEMV BLASFUNC(xgemv)
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM "
#تعریف_کرو GEMV BLASFUNC(zgemv)
#وإلا
#تعریف_کرو ERROR_NAME "CGEMM "
#تعریف_کرو GEMV BLASFUNC(cgemv)
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM3M "
#وإلا
#تعریف_کرو ERROR_NAME "CGEMM3M "
#جے_انت
#جے_انت
#جے_انت

#جے_بنا_تعریف GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#جے_انت

ساكن صحيح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#إذا !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#جے_انت
#وإلا
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#جے_انت
#جے_انت
};

#إذا defined(SMALL_میٹرکس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#تعریف_کرو USE_SMALL_میٹرکس_OPT 1
#وإلا
#تعریف_کرو USE_SMALL_میٹرکس_OPT 0
#جے_انت

#إذا USE_SMALL_میٹرکس_OPT
#جے_بنا_تعریف DYNAMIC_ARCH
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(جدول[idx]))
#وإلا
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#جے_انت


#جے_بنا_تعریف COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کرو GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کرو GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#وإلا

ساكن size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ساكن size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#تعریف_کرو ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کرو ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#جے_انت
#جے_انت

#إذا defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#تعریف_کرو XFEATURE_XTILEDATA 18
#تعریف_کرو ARCH_REQ_XCOMP_PERM 0x1023
ساكن صحيح openblas_amxtile_permission = 0;
ساكن صحيح init_amxtile_permission() {
  لما حالت =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  إذا (حالت != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    أرجع -1;
  }
  openblas_amxtile_permission = 1;
  أرجع 0;
}
#جے_انت

#جے_تعریف_شدہ DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#جے_انت

#جے_تعریف_شدہ SMP
#إذا defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ساكن مضمن صحيح get_gemm_optimal_nthreads_neoversev1(مزدوج MNK, صحيح ncpu) {
  أرجع
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
#جے_انت

#إذا defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ساكن مضمن صحيح get_gemm_optimal_nthreads_neoversev2(مزدوج MNK, صحيح ncpu) {
  أرجع
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
#جے_انت

ساكن مضمن صحيح get_gemm_optimal_nthreads(مزدوج MNK) {
  صحيح ncpu = عدد_cpu_avail(3);
#إذا defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  أرجع get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#وإلاإذا defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  أرجع get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#وإلاإذا defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  إذا (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    أرجع get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  إذا (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    أرجع get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#جے_انت
  إذا ( MNK <= (SMP_حد_MIN  * (مزدوج) GEMM_MULTITHREAD_حد) ) {
    أرجع 1;
  }
  وإلا {
    إذا (MNK/ncpu < SMP_حد_MIN*(مزدوج)GEMM_MULTITHREAD_حد) {
      أرجع MNK/(SMP_حد_MIN*(مزدوج)GEMM_MULTITHREAD_حد);
    }
    وإلا {
      أرجع ncpu;
    }
  }
}
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(حرف *TRANSA, حرف *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t دلیلاں;

  صحيح transa, transb, nrowa, nrowb;
  blasint info;

  حرف transA, transB;
  IFLOAT *بفر;
  IFLOAT *sa, *sb;

#جے_تعریف_شدہ SMP
  مزدوج MNK;
#إذا defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت
#جے_انت

#إذا defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  صحيح nodes;
#جے_انت

  PRINT_DEBUG_NAME;

  دلیلاں.m = *M;
  دلیلاں.n = *N;
  دلیلاں.k = *K;

  دلیلاں.a = (فراغ *)a;
  دلیلاں.b = (فراغ *)b;
  دلیلاں.c = (فراغ *)c;

  دلیلاں.lda = *ldA;
  دلیلاں.ldb = *ldB;
  دلیلاں.ldc = *ldC;

  دلیلاں.alpha = (فراغ *)alpha;
  دلیلاں.beta  = (فراغ *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  إذا (transA == 'N') transa = 0;
  إذا (transA == 'T') transa = 1;
#جے_بنا_تعریف COMPLEX
  إذا (transA == 'R') transa = 0;
  إذا (transA == 'C') transa = 1;
#وإلا
  إذا (transA == 'R') transa = 2;
  إذا (transA == 'C') transa = 3;
#جے_انت

  إذا (transB == 'N') transb = 0;
  إذا (transB == 'T') transb = 1;
#جے_بنا_تعریف COMPLEX
  إذا (transB == 'R') transb = 0;
  إذا (transB == 'C') transb = 1;
#وإلا
  إذا (transB == 'R') transb = 2;
  إذا (transB == 'C') transb = 3;
#جے_انت
  nrowa = دلیلاں.m;
  إذا (transa & 1) nrowa = دلیلاں.k;
  nrowb = دلیلاں.k;
  إذا (transb & 1) nrowb = دلیلاں.n;

  info = 0;

  إذا (دلیلاں.ldc < دلیلاں.m) info = 13;
  إذا (دلیلاں.ldb < nrowb)  info = 10;
  إذا (دلیلاں.lda < nrowa)  info =  8;
  إذا (دلیلاں.k < 0)        info =  5;
  إذا (دلیلاں.n < 0)        info =  4;
  إذا (دلیلاں.m < 0)        info =  3;
  إذا (transb < 0)        info =  2;
  إذا (transa < 0)        info =  1;

  إذا (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#جے_بنا_تعریف COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
	   فراغ *vbeta,
	   فراغ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#جے_انت

  blas_arg_t دلیلاں;
  صحيح transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#جے_تعریف_شدہ SMP
  مزدوج MNK;
#إذا defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت
#جے_انت

#إذا defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  صحيح nodes;
#جے_انت

  PRINT_DEBUG_CNAME;

#جے_بنا_تعریف COMPLEX
  دلیلاں.alpha = (فراغ *)&alpha;
  دلیلاں.beta  = (فراغ *)&beta;
#وإلا
  دلیلاں.alpha = (فراغ *)alpha;
  دلیلاں.beta  = (فراغ *)beta;
#جے_انت

  transa = -1;
  transb = -1;
  info   =  0;

  إذا (ترتیب == Cblas_خانہ_Major) {
    دلیلاں.m = m;
    دلیلاں.n = n;
    دلیلاں.k = k;

    دلیلاں.a = (فراغ *)a;
    دلیلاں.b = (فراغ *)b;
    دلیلاں.c = (فراغ *)c;

    دلیلاں.lda = lda;
    دلیلاں.ldb = ldb;
    دلیلاں.ldc = ldc;

    إذا (TransA == CblasNoTrans)     transa = 0;
    إذا (TransA == CblasTrans)       transa = 1;
#جے_بنا_تعریف COMPLEX
    إذا (TransA == CblasConjNoTrans) transa = 0;
    إذا (TransA == CblasConjTrans)   transa = 1;
#وإلا
    إذا (TransA == CblasConjNoTrans) transa = 2;
    إذا (TransA == CblasConjTrans)   transa = 3;
#جے_انت
    إذا (TransB == CblasNoTrans)     transb = 0;
    إذا (TransB == CblasTrans)       transb = 1;
#جے_بنا_تعریف COMPLEX
    إذا (TransB == CblasConjNoTrans) transb = 0;
    إذا (TransB == CblasConjTrans)   transb = 1;
#وإلا
    إذا (TransB == CblasConjNoTrans) transb = 2;
    إذا (TransB == CblasConjTrans)   transb = 3;
#جے_انت

    nrowa = دلیلاں.m;
    إذا (transa & 1) nrowa = دلیلاں.k;
    nrowb = دلیلاں.k;
    إذا (transb & 1) nrowb = دلیلاں.n;

    info = -1;

    إذا (دلیلاں.ldc < دلیلاں.m) info = 13;
    إذا (دلیلاں.ldb < nrowb)  info = 10;
    إذا (دلیلاں.lda < nrowa)  info =  8;
    إذا (دلیلاں.k < 0)        info =  5;
    إذا (دلیلاں.n < 0)        info =  4;
    إذا (دلیلاں.m < 0)        info =  3;
    إذا (transb < 0)        info =  2;
    إذا (transa < 0)        info =  1;
  }

  إذا (ترتیب == Cblas_صف_Major) {
    دلیلاں.m = n;
    دلیلاں.n = m;
    دلیلاں.k = k;

    دلیلاں.a = (فراغ *)b;
    دلیلاں.b = (فراغ *)a;
    دلیلاں.c = (فراغ *)c;

    دلیلاں.lda = ldb;
    دلیلاں.ldb = lda;
    دلیلاں.ldc = ldc;

    إذا (TransB == CblasNoTrans)     transa = 0;
    إذا (TransB == CblasTrans)       transa = 1;
#جے_بنا_تعریف COMPLEX
    إذا (TransB == CblasConjNoTrans) transa = 0;
    إذا (TransB == CblasConjTrans)   transa = 1;
#وإلا
    إذا (TransB == CblasConjNoTrans) transa = 2;
    إذا (TransB == CblasConjTrans)   transa = 3;
#جے_انت
    إذا (TransA == CblasNoTrans)     transb = 0;
    إذا (TransA == CblasTrans)       transb = 1;
#جے_بنا_تعریف COMPLEX
    إذا (TransA == CblasConjNoTrans) transb = 0;
    إذا (TransA == CblasConjTrans)   transb = 1;
#وإلا
    إذا (TransA == CblasConjNoTrans) transb = 2;
    إذا (TransA == CblasConjTrans)   transb = 3;
#جے_انت

    nrowa = دلیلاں.m;
    إذا (transa & 1) nrowa = دلیلاں.k;
    nrowb = دلیلاں.k;
    إذا (transb & 1) nrowb = دلیلاں.n;

    info = -1;

    إذا (دلیلاں.ldc < دلیلاں.m) info = 13;
    إذا (دلیلاں.ldb < nrowb)  info = 10;
    إذا (دلیلاں.lda < nrowa)  info =  8;
    إذا (دلیلاں.k < 0)        info =  5;
    إذا (دلیلاں.n < 0)        info =  4;
    إذا (دلیلاں.m < 0)        info =  3;
    إذا (transb < 0)        info =  2;
    إذا (transa < 0)        info =  1;

  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }


  إذا ((دلیلاں.m == 0) || (دلیلاں.n == 0)) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_سدھا)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
  إذا (support_avx512() )
#جے_انت
  إذا (ترتیب == Cblas_صف_Major && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_سدھا_PERFORMANT(m,n,k)) {
        SGEMM_سدھا(m, n, k, a, lda, b, ldb, c, ldc);
        أرجع;
  }
#جے_انت
#إذا defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_سدھا)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جے_انت
)
#جے_انت

  إذا (ترتیب == Cblas_صف_Major && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_سدھا_PERFORMANT(m,n,k)) {
        SGEMM_سدھا(m, n, k, a, lda, b, ldb, c, ldc);
        أرجع;
  }
وإلا إذا (ترتیب == Cblas_صف_Major && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_سدھا_PERFORMANT(m,n,k)) {
        SGEMM_سدھا_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
  }

#جے_انت //defined arm64
#جے_انت //defined complex


#جے_انت //defined CBLAS

#إذا !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_سدھا)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جے_انت
)
#جے_انت //defined dynarch
{
حرف* TA,*TB;
      إذا (transa & 1) 
        TA = "T";
	وإلا 
	TA= "N";
      إذا (transb & 1) 
        TB = "T";
	وإلا 
	TB= "N";
#جے_بنا_تعریف COMPLEX
      إذا (transa == 3) 
	TA= "T";
      إذا (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)دلیلاں.alpha;
FLOAT* be=(FLOAT*)دلیلاں.beta;
#جے_بنا_تعریف DOUBLE
        SME_SGEMM_KERNEL(TA,TB, دلیلاں.m, دلیلاں.n, دلیلاں.k, al, دلیلاں.a, دلیلاں.lda, دلیلاں.b, دلیلاں.ldb, be, دلیلاں.c, دلیلاں.ldc);
#وإلا
        SME_DGEMM_KERNEL(TA,TB, دلیلاں.m, دلیلاں.n, دلیلاں.k, al, دلیلاں.a, دلیلاں.lda, دلیلاں.b, دلیلاں.ldb, be, دلیلاں.c, دلیلاں.ldc);
#جے_انت
#وإلا
      إذا (transa == 2) 
	TA= "R";
      إذا (transb == 2) 
	TB= "R";
      إذا (transa == 3) 
	TA= "C";
      إذا (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)دلیلاں.alpha;
FLOAT* be=(FLOAT*)دلیلاں.beta;
#جے_بنا_تعریف DOUBLE
عائم _Complex c_al={al[0],al[1]};
عائم _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, دلیلاں.m, دلیلاں.n, دلیلاں.k, al[0],al[1], دلیلاں.a, دلیلاں.lda, دلیلاں.b, دلیلاں.ldb, be[0],be[1], دلیلاں.c, دلیلاں.ldc);
#وإلا
مزدوج _Complex c_al={al[0],al[1]};
مزدوج _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, دلیلاں.m, دلیلاں.n, دلیلاں.k, al[0],al[1], دلیلاں.a, دلیلاں.lda, دلیلاں.b, دلیلاں.ldb, be[0],be[1], دلیلاں.c, دلیلاں.ldc);
#جے_انت
#جے_انت
	أرجع;
}
#جے_انت //defined arm64

#جے_انت //defined b/hfloat16 

#إذا defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#إذا defined(DYNAMIC_ARCH)
  إذا (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    أرجع;
  }
#جے_انت
#إذا !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  إذا (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    أرجع;
  }
#جے_انت
#جے_انت  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  إذا ((دلیلاں.m == 0) || (دلیلاں.n == 0)) أرجع;

#إذا 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 دلیلاں.m, دلیلاں.n, دلیلاں.k, دلیلاں.lda, دلیلاں.ldb, دلیلاں.ldc);
#جے_انت

#إذا (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#تعریف_کرو BFLOAT16_GEMM_GEMV_FORWARD 1
#وإلا
#تعریف_کرو BFLOAT16_GEMM_GEMV_FORWARD 0
#جے_انت
#إذا (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#تعریف_کرو HFLOAT16_GEMM_GEMV_FORWARD 1
#وإلا
#تعریف_کرو HFLOAT16_GEMM_GEMV_FORWARD 0
#جے_انت

#إذا defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#إذا defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#وإلا
  bool have_tuned_gemv = true;
#جے_انت
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ساكن ثابت حرف gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  إذا (دلیلاں.k != 0) {
    إذا (دلیلاں.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = دلیلاں.m;
      blasint n = دلیلاں.k;
      blasint lda = دلیلاں.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      صحيح trans_idx = transa | ((transb & 2) << 1);
      حرف NT = gemv_trans[trans_idx];
      إذا (transa & 1) {
        m = دلیلاں.k;
        n = دلیلاں.m;
      }
      إذا (transb & 1) {
        inc_x = دلیلاں.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      إذا (is_efficient_gemv) {
        GEMV(&NT, &m, &n, دلیلاں.alpha, دلیلاں.a, &lda, دلیلاں.b, &inc_x, دلیلاں.beta, دلیلاں.c, &inc_y);
        أرجع;
      }
    }
    إذا (دلیلاں.m == 1) {
      blasint inc_x = دلیلاں.lda;
      blasint inc_y = دلیلاں.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = دلیلاں.k;
      blasint n = دلیلاں.n;
      blasint ldb = دلیلاں.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      صحيح trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      حرف NT = gemv_trans[trans_idx];
      إذا (transa & 1) {
        inc_x = 1;
      }
      إذا (transb & 1) {
        m = دلیلاں.n;
        n = دلیلاں.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      إذا (is_efficient_gemv) {
        GEMV(&NT, &m, &n, دلیلاں.alpha, دلیلاں.b, &ldb, دلیلاں.a, &inc_x, دلیلاں.beta, دلیلاں.c, &inc_y);
        أرجع;
      }
    }
  }
#جے_انت

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#إذا USE_SMALL_میٹرکس_OPT
إذا (*(FLOAT*)دلیلاں.alpha != 0.0) {
#إذا !defined(COMPLEX)
  إذا(GEMM_SMALL_میٹرکس_PERMIT(transa, transb, دلیلاں.m, دلیلاں.n, دلیلاں.k, *(FLOAT *)(دلیلاں.alpha), *(FLOAT *)(دلیلاں.beta))){

	  إذا(*(FLOAT *)(دلیلاں.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(دلیلاں.m, دلیلاں.n, دلیلاں.k, دلیلاں.a, دلیلاں.lda, *(FLOAT *)(دلیلاں.alpha), دلیلاں.b, دلیلاں.ldb, دلیلاں.c, دلیلاں.ldc);
	  }وإلا{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(دلیلاں.m, دلیلاں.n, دلیلاں.k, دلیلاں.a, دلیلاں.lda, *(FLOAT *)(دلیلاں.alpha), دلیلاں.b, دلیلاں.ldb, *(FLOAT *)(دلیلاں.beta), دلیلاں.c, دلیلاں.ldc);
	  }
	  أرجع;
  }
#وإلا
  إذا(GEMM_SMALL_میٹرکس_PERMIT(transa, transb, دلیلاں.m, دلیلاں.n, دلیلاں.k, alpha[0], alpha[1], beta[0], beta[1])){
	  إذا(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(دلیلاں.m, دلیلاں.n, دلیلاں.k, دلیلاں.a, دلیلاں.lda, alpha[0], alpha[1], دلیلاں.b, دلیلاں.ldb, دلیلاں.c, دلیلاں.ldc);
	  }وإلا{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(دلیلاں.m, دلیلاں.n, دلیلاں.k, دلیلاں.a, دلیلاں.lda, alpha[0], alpha[1], دلیلاں.b, دلیلاں.ldb, beta[0], beta[1], دلیلاں.c, دلیلاں.ldc);
	  }
	  أرجع;
  }
#جے_انت
}
#جے_انت

  بفر = (XFLOAT *)blas_حافظہ_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#إذا defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)بفر + (WhereAmI() & 0xf) * GEMM_وتھ_A);
#وإلا
  sa = (XFLOAT *)((BLASLONG)بفر +GEMM_وتھ_A);
#جے_انت
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ناپ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وتھ_B);

#جے_تعریف_شدہ SMP
#إذا defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  منوال |= (transa << BLAS_TRANSA_کھسکاؤ);
  منوال |= (transb << BLAS_TRANSB_کھسکاؤ);
#جے_انت

  MNK = (مزدوج) دلیلاں.m * (مزدوج) دلیلاں.n * (مزدوج) دلیلاں.k;
  دلیلاں.nthreads = get_gemm_optimal_nthreads(MNK);

  دلیلاں.common = عدم;

 إذا (دلیلاں.nthreads == 1) {
#جے_انت

    (gemm[(transb << 2) | transa])(&دلیلاں, عدم, عدم, sa, sb, 0);

#جے_تعریف_شدہ SMP

  } وإلا {

#جے_بنا_تعریف USE_SIMPLE_THREADED_LEVEL3

#جے_بنا_تعریف NO_AFFINITY
      nodes = get_عدد_nodes();

      إذا ((nodes > 1) && get_گرہ_equal()) {

	دلیلاں.nthreads /= nodes;

	gemm_thread_mn(منوال, &دلیلاں, عدم, عدم, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } وإلا {
#جے_انت

	(gemm[16 | (transb << 2) | transa])(&دلیلاں, عدم, عدم, sa, sb, 0);

#وإلا

	GEMM_THREAD(منوال, &دلیلاں, عدم, عدم, gemm[(transb << 2) | transa], sa, sb, دلیلاں.nthreads);

#جے_انت

#جے_بنا_تعریف USE_SIMPLE_THREADED_LEVEL3
#جے_بنا_تعریف NO_AFFINITY
      }
#جے_انت
#جے_انت

#جے_انت

#جے_تعریف_شدہ SMP
  }
#جے_انت

 blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_سرا(COMPSIZE * COMPSIZE, دلیلاں.m * دلیلاں.k + دلیلاں.k * دلیلاں.n + دلیلاں.m * دلیلاں.n, 2 * دلیلاں.m * دلیلاں.n * دلیلاں.k);

  IDEBUG_سرا;

  أرجع;
}
