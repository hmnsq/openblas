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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_نه_تعریف_شوی COMPLEX
#تعریف_کړه SMP_پوله_MIN 65536.0
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QGEMM "
#تعریف_کړه GEMV BLASFUNC(qgemv)
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DGEMM "
#تعریف_کړه GEMV BLASFUNC(dgemv)
#وإلاإذا defined(BFLOAT16)
#که_تعریف_شوی BGEMM
#تعریف_کړه ERROR_NAME "BGEMM "
#تعریف_کړه GEMV BLASFUNC(bgemv)
#وإلا
#تعریف_کړه ERROR_NAME "SBGEMM "
#تعریف_کړه GEMV BLASFUNC(sbgemv)
#د_که_پای
#وإلاإذا defined(HFLOAT16)
#تعریف_کړه ERROR_NAME "SHGEMM "
#وإلا
#تعریف_کړه ERROR_NAME "SGEMM "
#تعریف_کړه GEMV BLASFUNC(sgemv)
#د_که_پای
#وإلا
#تعریف_کړه SMP_پوله_MIN 8192.0
#که_نه_تعریف_شوی GEMM3M
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XGEMM "
#تعریف_کړه GEMV BLASFUNC(xgemv)
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZGEMM "
#تعریف_کړه GEMV BLASFUNC(zgemv)
#وإلا
#تعریف_کړه ERROR_NAME "CGEMM "
#تعریف_کړه GEMV BLASFUNC(cgemv)
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XGEMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZGEMM3M "
#وإلا
#تعریف_کړه ERROR_NAME "CGEMM3M "
#د_که_پای
#د_که_پای
#د_که_پای

#که_نه_تعریف_شوی GEMM_MULTITHREAD_پوله
#تعریف_کړه GEMM_MULTITHREAD_پوله 4
#د_که_پای

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
#د_که_پای
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
#د_که_پای
#د_که_پای
};

#إذا defined(SMALL_ماتریکس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#تعریف_کړه USE_SMALL_ماتریکس_OPT 1
#وإلا
#تعریف_کړه USE_SMALL_ماتریکس_OPT 0
#د_که_پای

#إذا USE_SMALL_ماتریکس_OPT
#که_نه_تعریف_شوی DYNAMIC_ARCH
#تعریف_کړه SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(جدول[idx]))
#وإلا
#تعریف_کړه SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#د_که_پای


#که_نه_تعریف_شوی COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کړه GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کړه GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
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

#تعریف_کړه ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کړه ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#د_که_پای
#د_که_پای

#إذا defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#تعریف_کړه XFEATURE_XTILEDATA 18
#تعریف_کړه ARCH_REQ_XCOMP_PERM 0x1023
ساكن صحيح openblas_amxtile_permission = 0;
ساكن صحيح init_amxtile_permission() {
  اوږد حالت =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  إذا (حالت != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    أرجع -1;
  }
  openblas_amxtile_permission = 1;
  أرجع 0;
}
#د_که_پای

#که_تعریف_شوی DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#د_که_پای

#که_تعریف_شوی SMP
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
#د_که_پای

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
#د_که_پای

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
#د_که_پای
  إذا ( MNK <= (SMP_پوله_MIN  * (مزدوج) GEMM_MULTITHREAD_پوله) ) {
    أرجع 1;
  }
  وإلا {
    إذا (MNK/ncpu < SMP_پوله_MIN*(مزدوج)GEMM_MULTITHREAD_پوله) {
      أرجع MNK/(SMP_پوله_MIN*(مزدوج)GEMM_MULTITHREAD_پوله);
    }
    وإلا {
      أرجع ncpu;
    }
  }
}
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *TRANSA, حرف *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t دلیلونه;

  صحيح transa, transb, nrowa, nrowb;
  blasint info;

  حرف transA, transB;
  IFLOAT *بفر;
  IFLOAT *sa, *sb;

#که_تعریف_شوی SMP
  مزدوج MNK;
#إذا defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای
#د_که_پای

#إذا defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  صحيح nodes;
#د_که_پای

  PRINT_DEBUG_NAME;

  دلیلونه.m = *M;
  دلیلونه.n = *N;
  دلیلونه.k = *K;

  دلیلونه.a = (فراغ *)a;
  دلیلونه.b = (فراغ *)b;
  دلیلونه.c = (فراغ *)c;

  دلیلونه.lda = *ldA;
  دلیلونه.ldb = *ldB;
  دلیلونه.ldc = *ldC;

  دلیلونه.alpha = (فراغ *)alpha;
  دلیلونه.beta  = (فراغ *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  إذا (transA == 'N') transa = 0;
  إذا (transA == 'T') transa = 1;
#که_نه_تعریف_شوی COMPLEX
  إذا (transA == 'R') transa = 0;
  إذا (transA == 'C') transa = 1;
#وإلا
  إذا (transA == 'R') transa = 2;
  إذا (transA == 'C') transa = 3;
#د_که_پای

  إذا (transB == 'N') transb = 0;
  إذا (transB == 'T') transb = 1;
#که_نه_تعریف_شوی COMPLEX
  إذا (transB == 'R') transb = 0;
  إذا (transB == 'C') transb = 1;
#وإلا
  إذا (transB == 'R') transb = 2;
  إذا (transB == 'C') transb = 3;
#د_که_پای
  nrowa = دلیلونه.m;
  إذا (transa & 1) nrowa = دلیلونه.k;
  nrowb = دلیلونه.k;
  إذا (transb & 1) nrowb = دلیلونه.n;

  info = 0;

  إذا (دلیلونه.ldc < دلیلونه.m) info = 13;
  إذا (دلیلونه.ldb < nrowb)  info = 10;
  إذا (دلیلونه.lda < nrowa)  info =  8;
  إذا (دلیلونه.k < 0)        info =  5;
  إذا (دلیلونه.n < 0)        info =  4;
  إذا (دلیلونه.m < 0)        info =  3;
  إذا (transb < 0)        info =  2;
  إذا (transa < 0)        info =  1;

  إذا (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#که_نه_تعریف_شوی COMPLEX
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
#د_که_پای

  blas_arg_t دلیلونه;
  صحيح transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#که_تعریف_شوی SMP
  مزدوج MNK;
#إذا defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای
#د_که_پای

#إذا defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  صحيح nodes;
#د_که_پای

  PRINT_DEBUG_CNAME;

#که_نه_تعریف_شوی COMPLEX
  دلیلونه.alpha = (فراغ *)&alpha;
  دلیلونه.beta  = (فراغ *)&beta;
#وإلا
  دلیلونه.alpha = (فراغ *)alpha;
  دلیلونه.beta  = (فراغ *)beta;
#د_که_پای

  transa = -1;
  transb = -1;
  info   =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    دلیلونه.m = m;
    دلیلونه.n = n;
    دلیلونه.k = k;

    دلیلونه.a = (فراغ *)a;
    دلیلونه.b = (فراغ *)b;
    دلیلونه.c = (فراغ *)c;

    دلیلونه.lda = lda;
    دلیلونه.ldb = ldb;
    دلیلونه.ldc = ldc;

    إذا (TransA == CblasNoTrans)     transa = 0;
    إذا (TransA == CblasTrans)       transa = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (TransA == CblasConjNoTrans) transa = 0;
    إذا (TransA == CblasConjTrans)   transa = 1;
#وإلا
    إذا (TransA == CblasConjNoTrans) transa = 2;
    إذا (TransA == CblasConjTrans)   transa = 3;
#د_که_پای
    إذا (TransB == CblasNoTrans)     transb = 0;
    إذا (TransB == CblasTrans)       transb = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (TransB == CblasConjNoTrans) transb = 0;
    إذا (TransB == CblasConjTrans)   transb = 1;
#وإلا
    إذا (TransB == CblasConjNoTrans) transb = 2;
    إذا (TransB == CblasConjTrans)   transb = 3;
#د_که_پای

    nrowa = دلیلونه.m;
    إذا (transa & 1) nrowa = دلیلونه.k;
    nrowb = دلیلونه.k;
    إذا (transb & 1) nrowb = دلیلونه.n;

    info = -1;

    إذا (دلیلونه.ldc < دلیلونه.m) info = 13;
    إذا (دلیلونه.ldb < nrowb)  info = 10;
    إذا (دلیلونه.lda < nrowa)  info =  8;
    إذا (دلیلونه.k < 0)        info =  5;
    إذا (دلیلونه.n < 0)        info =  4;
    إذا (دلیلونه.m < 0)        info =  3;
    إذا (transb < 0)        info =  2;
    إذا (transa < 0)        info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    دلیلونه.m = n;
    دلیلونه.n = m;
    دلیلونه.k = k;

    دلیلونه.a = (فراغ *)b;
    دلیلونه.b = (فراغ *)a;
    دلیلونه.c = (فراغ *)c;

    دلیلونه.lda = ldb;
    دلیلونه.ldb = lda;
    دلیلونه.ldc = ldc;

    إذا (TransB == CblasNoTrans)     transa = 0;
    إذا (TransB == CblasTrans)       transa = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (TransB == CblasConjNoTrans) transa = 0;
    إذا (TransB == CblasConjTrans)   transa = 1;
#وإلا
    إذا (TransB == CblasConjNoTrans) transa = 2;
    إذا (TransB == CblasConjTrans)   transa = 3;
#د_که_پای
    إذا (TransA == CblasNoTrans)     transb = 0;
    إذا (TransA == CblasTrans)       transb = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (TransA == CblasConjNoTrans) transb = 0;
    إذا (TransA == CblasConjTrans)   transb = 1;
#وإلا
    إذا (TransA == CblasConjNoTrans) transb = 2;
    إذا (TransA == CblasConjTrans)   transb = 3;
#د_که_پای

    nrowa = دلیلونه.m;
    إذا (transa & 1) nrowa = دلیلونه.k;
    nrowb = دلیلونه.k;
    إذا (transb & 1) nrowb = دلیلونه.n;

    info = -1;

    إذا (دلیلونه.ldc < دلیلونه.m) info = 13;
    إذا (دلیلونه.ldb < nrowb)  info = 10;
    إذا (دلیلونه.lda < nrowa)  info =  8;
    إذا (دلیلونه.k < 0)        info =  5;
    إذا (دلیلونه.n < 0)        info =  4;
    إذا (دلیلونه.m < 0)        info =  3;
    إذا (transb < 0)        info =  2;
    إذا (transa < 0)        info =  1;

  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }


  إذا ((دلیلونه.m == 0) || (دلیلونه.n == 0)) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
  إذا (support_avx512() )
#د_که_پای
  إذا (ترتیب == Cblas_لیکه_Major && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_مستقیم_PERFORMANT(m,n,k)) {
        SGEMM_مستقیم(m, n, k, a, lda, b, ldb, c, ldc);
        أرجع;
  }
#د_که_پای
#إذا defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#د_که_پای
)
#د_که_پای

  إذا (ترتیب == Cblas_لیکه_Major && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_مستقیم_PERFORMANT(m,n,k)) {
        SGEMM_مستقیم(m, n, k, a, lda, b, ldb, c, ldc);
        أرجع;
  }
وإلا إذا (ترتیب == Cblas_لیکه_Major && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_مستقیم_PERFORMANT(m,n,k)) {
        SGEMM_مستقیم_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
  }

#د_که_پای //defined arm64
#د_که_پای //defined complex


#د_که_پای //defined CBLAS

#إذا !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#د_که_پای
)
#د_که_پای //defined dynarch
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
#که_نه_تعریف_شوی COMPLEX
      إذا (transa == 3) 
	TA= "T";
      إذا (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)دلیلونه.alpha;
FLOAT* be=(FLOAT*)دلیلونه.beta;
#که_نه_تعریف_شوی DOUBLE
        SME_SGEMM_KERNEL(TA,TB, دلیلونه.m, دلیلونه.n, دلیلونه.k, al, دلیلونه.a, دلیلونه.lda, دلیلونه.b, دلیلونه.ldb, be, دلیلونه.c, دلیلونه.ldc);
#وإلا
        SME_DGEMM_KERNEL(TA,TB, دلیلونه.m, دلیلونه.n, دلیلونه.k, al, دلیلونه.a, دلیلونه.lda, دلیلونه.b, دلیلونه.ldb, be, دلیلونه.c, دلیلونه.ldc);
#د_که_پای
#وإلا
      إذا (transa == 2) 
	TA= "R";
      إذا (transb == 2) 
	TB= "R";
      إذا (transa == 3) 
	TA= "C";
      إذا (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)دلیلونه.alpha;
FLOAT* be=(FLOAT*)دلیلونه.beta;
#که_نه_تعریف_شوی DOUBLE
عائم _Complex c_al={al[0],al[1]};
عائم _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, دلیلونه.m, دلیلونه.n, دلیلونه.k, al[0],al[1], دلیلونه.a, دلیلونه.lda, دلیلونه.b, دلیلونه.ldb, be[0],be[1], دلیلونه.c, دلیلونه.ldc);
#وإلا
مزدوج _Complex c_al={al[0],al[1]};
مزدوج _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, دلیلونه.m, دلیلونه.n, دلیلونه.k, al[0],al[1], دلیلونه.a, دلیلونه.lda, دلیلونه.b, دلیلونه.ldb, be[0],be[1], دلیلونه.c, دلیلونه.ldc);
#د_که_پای
#د_که_پای
	أرجع;
}
#د_که_پای //defined arm64

#د_که_پای //defined b/hfloat16 

#إذا defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#إذا defined(DYNAMIC_ARCH)
  إذا (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    أرجع;
  }
#د_که_پای
#إذا !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  إذا (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    أرجع;
  }
#د_که_پای
#د_که_پای  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  إذا ((دلیلونه.m == 0) || (دلیلونه.n == 0)) أرجع;

#إذا 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 دلیلونه.m, دلیلونه.n, دلیلونه.k, دلیلونه.lda, دلیلونه.ldb, دلیلونه.ldc);
#د_که_پای

#إذا (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#تعریف_کړه BFLOAT16_GEMM_GEMV_FORWARD 1
#وإلا
#تعریف_کړه BFLOAT16_GEMM_GEMV_FORWARD 0
#د_که_پای
#إذا (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#تعریف_کړه HFLOAT16_GEMM_GEMV_FORWARD 1
#وإلا
#تعریف_کړه HFLOAT16_GEMM_GEMV_FORWARD 0
#د_که_پای

#إذا defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#إذا defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#وإلا
  bool have_tuned_gemv = true;
#د_که_پای
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ساكن ثابت حرف gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  إذا (دلیلونه.k != 0) {
    إذا (دلیلونه.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = دلیلونه.m;
      blasint n = دلیلونه.k;
      blasint lda = دلیلونه.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      صحيح trans_idx = transa | ((transb & 2) << 1);
      حرف NT = gemv_trans[trans_idx];
      إذا (transa & 1) {
        m = دلیلونه.k;
        n = دلیلونه.m;
      }
      إذا (transb & 1) {
        inc_x = دلیلونه.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      إذا (is_efficient_gemv) {
        GEMV(&NT, &m, &n, دلیلونه.alpha, دلیلونه.a, &lda, دلیلونه.b, &inc_x, دلیلونه.beta, دلیلونه.c, &inc_y);
        أرجع;
      }
    }
    إذا (دلیلونه.m == 1) {
      blasint inc_x = دلیلونه.lda;
      blasint inc_y = دلیلونه.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = دلیلونه.k;
      blasint n = دلیلونه.n;
      blasint ldb = دلیلونه.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      صحيح trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      حرف NT = gemv_trans[trans_idx];
      إذا (transa & 1) {
        inc_x = 1;
      }
      إذا (transb & 1) {
        m = دلیلونه.n;
        n = دلیلونه.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      إذا (is_efficient_gemv) {
        GEMV(&NT, &m, &n, دلیلونه.alpha, دلیلونه.b, &ldb, دلیلونه.a, &inc_x, دلیلونه.beta, دلیلونه.c, &inc_y);
        أرجع;
      }
    }
  }
#د_که_پای

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

#إذا USE_SMALL_ماتریکس_OPT
إذا (*(FLOAT*)دلیلونه.alpha != 0.0) {
#إذا !defined(COMPLEX)
  إذا(GEMM_SMALL_ماتریکس_PERMIT(transa, transb, دلیلونه.m, دلیلونه.n, دلیلونه.k, *(FLOAT *)(دلیلونه.alpha), *(FLOAT *)(دلیلونه.beta))){

	  إذا(*(FLOAT *)(دلیلونه.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(دلیلونه.m, دلیلونه.n, دلیلونه.k, دلیلونه.a, دلیلونه.lda, *(FLOAT *)(دلیلونه.alpha), دلیلونه.b, دلیلونه.ldb, دلیلونه.c, دلیلونه.ldc);
	  }وإلا{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(دلیلونه.m, دلیلونه.n, دلیلونه.k, دلیلونه.a, دلیلونه.lda, *(FLOAT *)(دلیلونه.alpha), دلیلونه.b, دلیلونه.ldb, *(FLOAT *)(دلیلونه.beta), دلیلونه.c, دلیلونه.ldc);
	  }
	  أرجع;
  }
#وإلا
  إذا(GEMM_SMALL_ماتریکس_PERMIT(transa, transb, دلیلونه.m, دلیلونه.n, دلیلونه.k, alpha[0], alpha[1], beta[0], beta[1])){
	  إذا(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(دلیلونه.m, دلیلونه.n, دلیلونه.k, دلیلونه.a, دلیلونه.lda, alpha[0], alpha[1], دلیلونه.b, دلیلونه.ldb, دلیلونه.c, دلیلونه.ldc);
	  }وإلا{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(دلیلونه.m, دلیلونه.n, دلیلونه.k, دلیلونه.a, دلیلونه.lda, alpha[0], alpha[1], دلیلونه.b, دلیلونه.ldb, beta[0], beta[1], دلیلونه.c, دلیلونه.ldc);
	  }
	  أرجع;
  }
#د_که_پای
}
#د_که_پای

  بفر = (XFLOAT *)blas_حافظه_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#إذا defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)بفر + (WhereAmI() & 0xf) * GEMM_دننه_واټن_A);
#وإلا
  sa = (XFLOAT *)((BLASLONG)بفر +GEMM_دننه_واټن_A);
#د_که_پای
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_دننه_واټن_B);

#که_تعریف_شوی SMP
#إذا defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  ډېرپېښ |= (transa << BLAS_TRANSA_ښويدنه);
  ډېرپېښ |= (transb << BLAS_TRANSB_ښويدنه);
#د_که_پای

  MNK = (مزدوج) دلیلونه.m * (مزدوج) دلیلونه.n * (مزدوج) دلیلونه.k;
  دلیلونه.nthreads = get_gemm_optimal_nthreads(MNK);

  دلیلونه.common = عدم;

 إذا (دلیلونه.nthreads == 1) {
#د_که_پای

    (gemm[(transb << 2) | transa])(&دلیلونه, عدم, عدم, sa, sb, 0);

#که_تعریف_شوی SMP

  } وإلا {

#که_نه_تعریف_شوی USE_SIMPLE_THREADED_LEVEL3

#که_نه_تعریف_شوی NO_AFFINITY
      nodes = get_عدد_nodes();

      إذا ((nodes > 1) && get_غوټه_equal()) {

	دلیلونه.nthreads /= nodes;

	gemm_thread_mn(ډېرپېښ, &دلیلونه, عدم, عدم, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } وإلا {
#د_که_پای

	(gemm[16 | (transb << 2) | transa])(&دلیلونه, عدم, عدم, sa, sb, 0);

#وإلا

	GEMM_THREAD(ډېرپېښ, &دلیلونه, عدم, عدم, gemm[(transb << 2) | transa], sa, sb, دلیلونه.nthreads);

#د_که_پای

#که_نه_تعریف_شوی USE_SIMPLE_THREADED_LEVEL3
#که_نه_تعریف_شوی NO_AFFINITY
      }
#د_که_پای
#د_که_پای

#د_که_پای

#که_تعریف_شوی SMP
  }
#د_که_پای

 blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE * COMPSIZE, دلیلونه.m * دلیلونه.k + دلیلونه.k * دلیلونه.n + دلیلونه.m * دلیلونه.n, 2 * دلیلونه.m * دلیلونه.n * دلیلونه.k);

  IDEBUG_پای_ټکی;

  أرجع;
}
