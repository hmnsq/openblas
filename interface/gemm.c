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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdlib.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdbool.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 65536.0
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(qgemv)
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(dgemv)
#ਜਾਂਜੇ defined(BFLOAT16)
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ BGEMM
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "BGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(bgemv)
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SBGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(sbgemv)
#ਜੇ_ਅੰਤ
#ਜਾਂਜੇ defined(HFLOAT16)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SHGEMM "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(sgemv)
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 8192.0
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ GEMM3M
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(xgemv)
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(zgemv)
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGEMM "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV BLASFUNC(cgemv)
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGEMM3M "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGEMM3M "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGEMM3M "
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ GEMM_MULTITHREAD_ਹੱਦ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_MULTITHREAD_ਹੱਦ 4
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#ਜੇ !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#ਜੇ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#ਜੇ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
};

#ਜੇ defined(SMALL_ਮੈਟ੍ਰਿਕਸ_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT 1
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT 0
#ਜੇ_ਅੰਤ

#ਜੇ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ DYNAMIC_ARCH
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMALL_KERNEL_ADDR(ਸਾਰਣੀ, idx) ((ਖਾਲੀ *)(ਸਾਰਣੀ[idx]))
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMALL_KERNEL_ADDR(ਸਾਰਣੀ, idx) ((ਖਾਲੀ *)(*(uintptr_t *)((ਅੱਖਰ *)gotoblas + (size_t)(ਸਾਰਣੀ[idx]))))
#ਜੇ_ਅੰਤ


#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
ਸਥਿਰ size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ਸਥਿਰ size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_SMALL_KERNEL_B0(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_SMALL_KERNEL(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ਨਹੀਂਤਾਂ

ਸਥਿਰ size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ਸਥਿਰ size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ZGEMM_SMALL_KERNEL(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ZGEMM_SMALL_KERNEL_B0(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ XFEATURE_XTILEDATA 18
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ARCH_REQ_XCOMP_PERM 0x1023
ਸਥਿਰ ਪੂਰਨ openblas_amxtile_permission = 0;
ਸਥਿਰ ਪੂਰਨ init_amxtile_permission() {
  ਲੰਮਾ ਸਥਿਤੀ =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  ਜੇ (ਸਥਿਤੀ != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    ਵਾਪਸ -1;
  }
  openblas_amxtile_permission = 1;
  ਵਾਪਸ 0;
}
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DYNAMIC_ARCH
ਬਾਹਰੀ ਅੱਖਰ* gotoblas_corename(ਖਾਲੀ);
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ਸਥਿਰ ਸਤਰਅੰਦਰ ਪੂਰਨ get_gemm_optimal_nthreads_neoversev1(ਦੁੱਗਣਾ MNK, ਪੂਰਨ ncpu) {
  ਵਾਪਸ
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
#ਜੇ_ਅੰਤ

#ਜੇ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ਸਥਿਰ ਸਤਰਅੰਦਰ ਪੂਰਨ get_gemm_optimal_nthreads_neoversev2(ਦੁੱਗਣਾ MNK, ਪੂਰਨ ncpu) {
  ਵਾਪਸ
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
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਸਤਰਅੰਦਰ ਪੂਰਨ get_gemm_optimal_nthreads(ਦੁੱਗਣਾ MNK) {
  ਪੂਰਨ ncpu = ਅੰਕ_cpu_avail(3);
#ਜੇ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ਵਾਪਸ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#ਜਾਂਜੇ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ਵਾਪਸ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#ਜਾਂਜੇ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ਜੇ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ਵਾਪਸ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  ਜੇ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ਵਾਪਸ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#ਜੇ_ਅੰਤ
  ਜੇ ( MNK <= (SMP_ਹੱਦ_MIN  * (ਦੁੱਗਣਾ) GEMM_MULTITHREAD_ਹੱਦ) ) {
    ਵਾਪਸ 1;
  }
  ਨਹੀਂਤਾਂ {
    ਜੇ (MNK/ncpu < SMP_ਹੱਦ_MIN*(ਦੁੱਗਣਾ)GEMM_MULTITHREAD_ਹੱਦ) {
      ਵਾਪਸ MNK/(SMP_ਹੱਦ_MIN*(ਦੁੱਗਣਾ)GEMM_MULTITHREAD_ਹੱਦ);
    }
    ਨਹੀਂਤਾਂ {
      ਵਾਪਸ ncpu;
    }
  }
}
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *TRANSA, ਅੱਖਰ *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t ਦਲੀਲਾਂ;

  ਪੂਰਨ transa, transb, nrowa, nrowb;
  blasint info;

  ਅੱਖਰ transA, transB;
  IFLOAT *ਬਫਰ;
  IFLOAT *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਦੁੱਗਣਾ MNK;
#ਜੇ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ਪੂਰਨ nodes;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_NAME;

  ਦਲੀਲਾਂ.m = *M;
  ਦਲੀਲਾਂ.n = *N;
  ਦਲੀਲਾਂ.k = *K;

  ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.b = (ਖਾਲੀ *)b;
  ਦਲੀਲਾਂ.c = (ਖਾਲੀ *)c;

  ਦਲੀਲਾਂ.lda = *ldA;
  ਦਲੀਲਾਂ.ldb = *ldB;
  ਦਲੀਲਾਂ.ldc = *ldC;

  ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)alpha;
  ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  ਜੇ (transA == 'N') transa = 0;
  ਜੇ (transA == 'T') transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (transA == 'R') transa = 0;
  ਜੇ (transA == 'C') transa = 1;
#ਨਹੀਂਤਾਂ
  ਜੇ (transA == 'R') transa = 2;
  ਜੇ (transA == 'C') transa = 3;
#ਜੇ_ਅੰਤ

  ਜੇ (transB == 'N') transb = 0;
  ਜੇ (transB == 'T') transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (transB == 'R') transb = 0;
  ਜੇ (transB == 'C') transb = 1;
#ਨਹੀਂਤਾਂ
  ਜੇ (transB == 'R') transb = 2;
  ਜੇ (transB == 'C') transb = 3;
#ਜੇ_ਅੰਤ
  nrowa = ਦਲੀਲਾਂ.m;
  ਜੇ (transa & 1) nrowa = ਦਲੀਲਾਂ.k;
  nrowb = ਦਲੀਲਾਂ.k;
  ਜੇ (transb & 1) nrowb = ਦਲੀਲਾਂ.n;

  info = 0;

  ਜੇ (ਦਲੀਲਾਂ.ldc < ਦਲੀਲਾਂ.m) info = 13;
  ਜੇ (ਦਲੀਲਾਂ.ldb < nrowb)  info = 10;
  ਜੇ (ਦਲੀਲਾਂ.lda < nrowa)  info =  8;
  ਜੇ (ਦਲੀਲਾਂ.k < 0)        info =  5;
  ਜੇ (ਦਲੀਲਾਂ.n < 0)        info =  4;
  ਜੇ (ਦਲੀਲਾਂ.m < 0)        info =  3;
  ਜੇ (transb < 0)        info =  2;
  ਜੇ (transa < 0)        info =  1;

  ਜੇ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_TRANSPOSE TransA, ਗਿਣਤੀ CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *valpha,
	   ਖਾਲੀ *va, blasint lda,
	   ਖਾਲੀ *vb, blasint ldb,
	   ਖਾਲੀ *vbeta,
	   ਖਾਲੀ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#ਜੇ_ਅੰਤ

  blas_arg_t ਦਲੀਲਾਂ;
  ਪੂਰਨ transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *ਬਫਰ;
  XFLOAT *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਦੁੱਗਣਾ MNK;
#ਜੇ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ਪੂਰਨ nodes;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)&alpha;
  ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)&beta;
#ਨਹੀਂਤਾਂ
  ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)alpha;
  ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)beta;
#ਜੇ_ਅੰਤ

  transa = -1;
  transb = -1;
  info   =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਦਲੀਲਾਂ.m = m;
    ਦਲੀਲਾਂ.n = n;
    ਦਲੀਲਾਂ.k = k;

    ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
    ਦਲੀਲਾਂ.b = (ਖਾਲੀ *)b;
    ਦਲੀਲਾਂ.c = (ਖਾਲੀ *)c;

    ਦਲੀਲਾਂ.lda = lda;
    ਦਲੀਲਾਂ.ldb = ldb;
    ਦਲੀਲਾਂ.ldc = ldc;

    ਜੇ (TransA == CblasNoTrans)     transa = 0;
    ਜੇ (TransA == CblasTrans)       transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (TransA == CblasConjNoTrans) transa = 0;
    ਜੇ (TransA == CblasConjTrans)   transa = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (TransA == CblasConjNoTrans) transa = 2;
    ਜੇ (TransA == CblasConjTrans)   transa = 3;
#ਜੇ_ਅੰਤ
    ਜੇ (TransB == CblasNoTrans)     transb = 0;
    ਜੇ (TransB == CblasTrans)       transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (TransB == CblasConjNoTrans) transb = 0;
    ਜੇ (TransB == CblasConjTrans)   transb = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (TransB == CblasConjNoTrans) transb = 2;
    ਜੇ (TransB == CblasConjTrans)   transb = 3;
#ਜੇ_ਅੰਤ

    nrowa = ਦਲੀਲਾਂ.m;
    ਜੇ (transa & 1) nrowa = ਦਲੀਲਾਂ.k;
    nrowb = ਦਲੀਲਾਂ.k;
    ਜੇ (transb & 1) nrowb = ਦਲੀਲਾਂ.n;

    info = -1;

    ਜੇ (ਦਲੀਲਾਂ.ldc < ਦਲੀਲਾਂ.m) info = 13;
    ਜੇ (ਦਲੀਲਾਂ.ldb < nrowb)  info = 10;
    ਜੇ (ਦਲੀਲਾਂ.lda < nrowa)  info =  8;
    ਜੇ (ਦਲੀਲਾਂ.k < 0)        info =  5;
    ਜੇ (ਦਲੀਲਾਂ.n < 0)        info =  4;
    ਜੇ (ਦਲੀਲਾਂ.m < 0)        info =  3;
    ਜੇ (transb < 0)        info =  2;
    ਜੇ (transa < 0)        info =  1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਦਲੀਲਾਂ.m = n;
    ਦਲੀਲਾਂ.n = m;
    ਦਲੀਲਾਂ.k = k;

    ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)b;
    ਦਲੀਲਾਂ.b = (ਖਾਲੀ *)a;
    ਦਲੀਲਾਂ.c = (ਖਾਲੀ *)c;

    ਦਲੀਲਾਂ.lda = ldb;
    ਦਲੀਲਾਂ.ldb = lda;
    ਦਲੀਲਾਂ.ldc = ldc;

    ਜੇ (TransB == CblasNoTrans)     transa = 0;
    ਜੇ (TransB == CblasTrans)       transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (TransB == CblasConjNoTrans) transa = 0;
    ਜੇ (TransB == CblasConjTrans)   transa = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (TransB == CblasConjNoTrans) transa = 2;
    ਜੇ (TransB == CblasConjTrans)   transa = 3;
#ਜੇ_ਅੰਤ
    ਜੇ (TransA == CblasNoTrans)     transb = 0;
    ਜੇ (TransA == CblasTrans)       transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (TransA == CblasConjNoTrans) transb = 0;
    ਜੇ (TransA == CblasConjTrans)   transb = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (TransA == CblasConjNoTrans) transb = 2;
    ਜੇ (TransA == CblasConjTrans)   transb = 3;
#ਜੇ_ਅੰਤ

    nrowa = ਦਲੀਲਾਂ.m;
    ਜੇ (transa & 1) nrowa = ਦਲੀਲਾਂ.k;
    nrowb = ਦਲੀਲਾਂ.k;
    ਜੇ (transb & 1) nrowb = ਦਲੀਲਾਂ.n;

    info = -1;

    ਜੇ (ਦਲੀਲਾਂ.ldc < ਦਲੀਲਾਂ.m) info = 13;
    ਜੇ (ਦਲੀਲਾਂ.ldb < nrowb)  info = 10;
    ਜੇ (ਦਲੀਲਾਂ.lda < nrowa)  info =  8;
    ਜੇ (ਦਲੀਲਾਂ.k < 0)        info =  5;
    ਜੇ (ਦਲੀਲਾਂ.n < 0)        info =  4;
    ਜੇ (ਦਲੀਲਾਂ.m < 0)        info =  3;
    ਜੇ (transb < 0)        info =  2;
    ਜੇ (transa < 0)        info =  1;

  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }


  ਜੇ ((ਦਲੀਲਾਂ.m == 0) || (ਦਲੀਲਾਂ.n == 0)) ਵਾਪਸ;
#ਜੇ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ਜੇ defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_ਸਿੱਧਾ)||defined(DYNAMIC_ARCH))
#ਜੇ defined(DYNAMIC_ARCH)
  ਜੇ (support_avx512() )
#ਜੇ_ਅੰਤ
  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ਸਿੱਧਾ_PERFORMANT(m,n,k)) {
        SGEMM_ਸਿੱਧਾ(m, n, k, a, lda, b, ldb, c, ldc);
        ਵਾਪਸ;
  }
#ਜੇ_ਅੰਤ
#ਜੇ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_ਸਿੱਧਾ)||defined(DYNAMIC_ARCH))
#ਜੇ defined(DYNAMIC_ARCH)
ਜੇ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ਜੇ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ਜੇ_ਅੰਤ
)
#ਜੇ_ਅੰਤ

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ਸਿੱਧਾ_PERFORMANT(m,n,k)) {
        SGEMM_ਸਿੱਧਾ(m, n, k, a, lda, b, ldb, c, ldc);
        ਵਾਪਸ;
  }
ਨਹੀਂਤਾਂ ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ਸਿੱਧਾ_PERFORMANT(m,n,k)) {
        SGEMM_ਸਿੱਧਾ_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ //defined arm64
#ਜੇ_ਅੰਤ //defined complex


#ਜੇ_ਅੰਤ //defined CBLAS

#ਜੇ !defined(BFLOAT16)  && !defined(HFLOAT16)
#ਜੇ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_ਸਿੱਧਾ)||defined(DYNAMIC_ARCH))
#ਜੇ defined(DYNAMIC_ARCH)
ਜੇ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ਜੇ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ਜੇ_ਅੰਤ
)
#ਜੇ_ਅੰਤ //defined dynarch
{
ਅੱਖਰ* TA,*TB;
      ਜੇ (transa & 1) 
        TA = "T";
	ਨਹੀਂਤਾਂ 
	TA= "N";
      ਜੇ (transb & 1) 
        TB = "T";
	ਨਹੀਂਤਾਂ 
	TB= "N";
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (transa == 3) 
	TA= "T";
      ਜੇ (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)ਦਲੀਲਾਂ.alpha;
FLOAT* be=(FLOAT*)ਦਲੀਲਾਂ.beta;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ DOUBLE
        SME_SGEMM_KERNEL(TA,TB, ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, al, ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, be, ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
#ਨਹੀਂਤਾਂ
        SME_DGEMM_KERNEL(TA,TB, ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, al, ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, be, ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
      ਜੇ (transa == 2) 
	TA= "R";
      ਜੇ (transb == 2) 
	TB= "R";
      ਜੇ (transa == 3) 
	TA= "C";
      ਜੇ (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)ਦਲੀਲਾਂ.alpha;
FLOAT* be=(FLOAT*)ਦਲੀਲਾਂ.beta;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ DOUBLE
ਦਸ਼ਮਲਵ _Complex c_al={al[0],al[1]};
ਦਸ਼ਮਲਵ _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, al[0],al[1], ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, be[0],be[1], ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
#ਨਹੀਂਤਾਂ
ਦੁੱਗਣਾ _Complex c_al={al[0],al[1]};
ਦੁੱਗਣਾ _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, al[0],al[1], ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, be[0],be[1], ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
	ਵਾਪਸ;
}
#ਜੇ_ਅੰਤ //defined arm64

#ਜੇ_ਅੰਤ //defined b/hfloat16 

#ਜੇ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ਜੇ defined(DYNAMIC_ARCH)
  ਜੇ (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ਵਾਪਸ;
  }
#ਜੇ_ਅੰਤ
#ਜੇ !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  ਜੇ (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ਵਾਪਸ;
  }
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  ਜੇ ((ਦਲੀਲਾਂ.m == 0) || (ਦਲੀਲਾਂ.n == 0)) ਵਾਪਸ;

#ਜੇ 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, ਦਲੀਲਾਂ.lda, ਦਲੀਲਾਂ.ldb, ਦਲੀਲਾਂ.ldc);
#ਜੇ_ਅੰਤ

#ਜੇ (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ BFLOAT16_GEMM_GEMV_FORWARD 1
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ BFLOAT16_GEMM_GEMV_FORWARD 0
#ਜੇ_ਅੰਤ
#ਜੇ (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ HFLOAT16_GEMM_GEMV_FORWARD 1
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ HFLOAT16_GEMM_GEMV_FORWARD 0
#ਜੇ_ਅੰਤ

#ਜੇ defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#ਜੇ defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#ਨਹੀਂਤਾਂ
  bool have_tuned_gemv = true;
#ਜੇ_ਅੰਤ
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ਸਥਿਰ ਅਚਲ ਅੱਖਰ gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  ਜੇ (ਦਲੀਲਾਂ.k != 0) {
    ਜੇ (ਦਲੀਲਾਂ.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = ਦਲੀਲਾਂ.m;
      blasint n = ਦਲੀਲਾਂ.k;
      blasint lda = ਦਲੀਲਾਂ.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      ਪੂਰਨ trans_idx = transa | ((transb & 2) << 1);
      ਅੱਖਰ NT = gemv_trans[trans_idx];
      ਜੇ (transa & 1) {
        m = ਦਲੀਲਾਂ.k;
        n = ਦਲੀਲਾਂ.m;
      }
      ਜੇ (transb & 1) {
        inc_x = ਦਲੀਲਾਂ.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      ਜੇ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, ਦਲੀਲਾਂ.alpha, ਦਲੀਲਾਂ.a, &lda, ਦਲੀਲਾਂ.b, &inc_x, ਦਲੀਲਾਂ.beta, ਦਲੀਲਾਂ.c, &inc_y);
        ਵਾਪਸ;
      }
    }
    ਜੇ (ਦਲੀਲਾਂ.m == 1) {
      blasint inc_x = ਦਲੀਲਾਂ.lda;
      blasint inc_y = ਦਲੀਲਾਂ.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = ਦਲੀਲਾਂ.k;
      blasint n = ਦਲੀਲਾਂ.n;
      blasint ldb = ਦਲੀਲਾਂ.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      ਪੂਰਨ trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      ਅੱਖਰ NT = gemv_trans[trans_idx];
      ਜੇ (transa & 1) {
        inc_x = 1;
      }
      ਜੇ (transb & 1) {
        m = ਦਲੀਲਾਂ.n;
        n = ਦਲੀਲਾਂ.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      ਜੇ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, ਦਲੀਲਾਂ.alpha, ਦਲੀਲਾਂ.b, &ldb, ਦਲੀਲਾਂ.a, &inc_x, ਦਲੀਲਾਂ.beta, ਦਲੀਲਾਂ.c, &inc_y);
        ਵਾਪਸ;
      }
    }
  }
#ਜੇ_ਅੰਤ

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

#ਜੇ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
ਜੇ (*(FLOAT*)ਦਲੀਲਾਂ.alpha != 0.0) {
#ਜੇ !defined(COMPLEX)
  ਜੇ(GEMM_SMALL_ਮੈਟ੍ਰਿਕਸ_PERMIT(transa, transb, ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, *(FLOAT *)(ਦਲੀਲਾਂ.alpha), *(FLOAT *)(ਦਲੀਲਾਂ.beta))){

	  ਜੇ(*(FLOAT *)(ਦਲੀਲਾਂ.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, *(FLOAT *)(ਦਲੀਲਾਂ.alpha), ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
	  }ਨਹੀਂਤਾਂ{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, *(FLOAT *)(ਦਲੀਲਾਂ.alpha), ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, *(FLOAT *)(ਦਲੀਲਾਂ.beta), ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
	  }
	  ਵਾਪਸ;
  }
#ਨਹੀਂਤਾਂ
  ਜੇ(GEMM_SMALL_ਮੈਟ੍ਰਿਕਸ_PERMIT(transa, transb, ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, alpha[0], alpha[1], beta[0], beta[1])){
	  ਜੇ(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, alpha[0], alpha[1], ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
	  }ਨਹੀਂਤਾਂ{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(ਦਲੀਲਾਂ.m, ਦਲੀਲਾਂ.n, ਦਲੀਲਾਂ.k, ਦਲੀਲਾਂ.a, ਦਲੀਲਾਂ.lda, alpha[0], alpha[1], ਦਲੀਲਾਂ.b, ਦਲੀਲਾਂ.ldb, beta[0], beta[1], ਦਲੀਲਾਂ.c, ਦਲੀਲਾਂ.ldc);
	  }
	  ਵਾਪਸ;
  }
#ਜੇ_ਅੰਤ
}
#ਜੇ_ਅੰਤ

  ਬਫਰ = (XFLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#ਜੇ defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)ਬਫਰ + (WhereAmI() & 0xf) * GEMM_ਵਿੱਥ_A);
#ਨਹੀਂਤਾਂ
  sa = (XFLOAT *)((BLASLONG)ਬਫਰ +GEMM_ਵਿੱਥ_A);
#ਜੇ_ਅੰਤ
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ਮਾਪ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ਵਿੱਥ_B);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  ਬਹੁਲਕ |= (transa << BLAS_TRANSA_ਖਿਸਕਾਓ);
  ਬਹੁਲਕ |= (transb << BLAS_TRANSB_ਖਿਸਕਾਓ);
#ਜੇ_ਅੰਤ

  MNK = (ਦੁੱਗਣਾ) ਦਲੀਲਾਂ.m * (ਦੁੱਗਣਾ) ਦਲੀਲਾਂ.n * (ਦੁੱਗਣਾ) ਦਲੀਲਾਂ.k;
  ਦਲੀਲਾਂ.nthreads = get_gemm_optimal_nthreads(MNK);

  ਦਲੀਲਾਂ.common = ਨਿਰਾ;

 ਜੇ (ਦਲੀਲਾਂ.nthreads == 1) {
#ਜੇ_ਅੰਤ

    (gemm[(transb << 2) | transa])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP

  } ਨਹੀਂਤਾਂ {

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ USE_SIMPLE_THREADED_LEVEL3

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ NO_AFFINITY
      nodes = get_ਅੰਕ_nodes();

      ਜੇ ((nodes > 1) && get_ਗੰਢ_equal()) {

	ਦਲੀਲਾਂ.nthreads /= nodes;

	gemm_thread_mn(ਬਹੁਲਕ, &ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } ਨਹੀਂਤਾਂ {
#ਜੇ_ਅੰਤ

	(gemm[16 | (transb << 2) | transa])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

#ਨਹੀਂਤਾਂ

	GEMM_THREAD(ਬਹੁਲਕ, &ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, gemm[(transb << 2) | transa], sa, sb, ਦਲੀਲਾਂ.nthreads);

#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ USE_SIMPLE_THREADED_LEVEL3
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ NO_AFFINITY
      }
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  }
#ਜੇ_ਅੰਤ

 blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE * COMPSIZE, ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.k + ਦਲੀਲਾਂ.k * ਦਲੀਲਾਂ.n + ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n, 2 * ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.k);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
