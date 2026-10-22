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
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 65536.0
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(qgemv)
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(dgemv)
#ಅಥವಾಆದರೆ defined(BFLOAT16)
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ BGEMM
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "BGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(bgemv)
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SBGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(sbgemv)
#ಆದರೆ_ಅಂತ್ಯ
#ಅಥವಾಆದರೆ defined(HFLOAT16)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SHGEMM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(sgemv)
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 8192.0
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM3M
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(xgemv)
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(zgemv)
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEMM "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV BLASFUNC(cgemv)
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGEMM3M "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEMM3M "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEMM3M "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM_MULTITHREAD_ಮಿತಿ
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_MULTITHREAD_ಮಿತಿ 4
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#ಆದರೆ !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#ಆದರೆ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
};

#ಆದರೆ defined(SMALL_ಮಾತೃಕೆ_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#ವ್ಯಾಖ್ಯಾನಿಸು USE_SMALL_ಮಾತೃಕೆ_OPT 1
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು USE_SMALL_ಮಾತೃಕೆ_OPT 0
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ USE_SMALL_ಮಾತೃಕೆ_OPT
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ DYNAMIC_ARCH
#ವ್ಯಾಖ್ಯಾನಿಸು SMALL_KERNEL_ADDR(ಕೋಷ್ಟಕ, idx) ((ಶೂನ್ಯ *)(ಕೋಷ್ಟಕ[idx]))
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMALL_KERNEL_ADDR(ಕೋಷ್ಟಕ, idx) ((ಶೂನ್ಯ *)(*(uintptr_t *)((ಅಕ್ಷರ *)gotoblas + (size_t)(ಕೋಷ್ಟಕ[idx]))))
#ಆದರೆ_ಅಂತ್ಯ


#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
ಸ್ಥಿರ size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ಸ್ಥಿರ size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ಇಲ್ಲದಿದ್ದರೆ

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

#ವ್ಯಾಖ್ಯಾನಿಸು ZGEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#ವ್ಯಾಖ್ಯಾನಿಸು ZGEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ವ್ಯಾಖ್ಯಾನಿಸು XFEATURE_XTILEDATA 18
#ವ್ಯಾಖ್ಯಾನಿಸು ARCH_REQ_XCOMP_PERM 0x1023
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ openblas_amxtile_permission = 0;
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ init_amxtile_permission() {
  ದೀರ್ಘ ಸ್ಥಿತಿ_2 =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  ಆದರೆ (ಸ್ಥಿತಿ_2 != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    ಹಿಂತಿರುಗಿ -1;
  }
  openblas_amxtile_permission = 1;
  ಹಿಂತಿರುಗಿ 0;
}
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DYNAMIC_ARCH
ಬಾಹ್ಯ ಅಕ್ಷರ* gotoblas_corename(ಶೂನ್ಯ);
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ಆದರೆ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ಸ್ಥಿರ ಸಾಲೊಳಗೆ ಪೂರ್ಣಾಂಕ get_gemm_optimal_nthreads_neoversev1(ದ್ವಿಗುಣ MNK, ಪೂರ್ಣಾಂಕ ncpu) {
  ಹಿಂತಿರುಗಿ
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
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ಸ್ಥಿರ ಸಾಲೊಳಗೆ ಪೂರ್ಣಾಂಕ get_gemm_optimal_nthreads_neoversev2(ದ್ವಿಗುಣ MNK, ಪೂರ್ಣಾಂಕ ncpu) {
  ಹಿಂತಿರುಗಿ
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
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಸಾಲೊಳಗೆ ಪೂರ್ಣಾಂಕ get_gemm_optimal_nthreads(ದ್ವಿಗುಣ MNK) {
  ಪೂರ್ಣಾಂಕ ncpu = ಅಂಕೆ_cpu_avail(3);
#ಆದರೆ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ಹಿಂತಿರುಗಿ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#ಅಥವಾಆದರೆ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ಹಿಂತಿರುಗಿ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#ಅಥವಾಆದರೆ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ಆದರೆ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ಹಿಂತಿರುಗಿ get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  ಆದರೆ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ಹಿಂತಿರುಗಿ get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#ಆದರೆ_ಅಂತ್ಯ
  ಆದರೆ ( MNK <= (SMP_ಮಿತಿ_MIN  * (ದ್ವಿಗುಣ) GEMM_MULTITHREAD_ಮಿತಿ) ) {
    ಹಿಂತಿರುಗಿ 1;
  }
  ಇಲ್ಲದಿದ್ದರೆ {
    ಆದರೆ (MNK/ncpu < SMP_ಮಿತಿ_MIN*(ದ್ವಿಗುಣ)GEMM_MULTITHREAD_ಮಿತಿ) {
      ಹಿಂತಿರುಗಿ MNK/(SMP_ಮಿತಿ_MIN*(ದ್ವಿಗುಣ)GEMM_MULTITHREAD_ಮಿತಿ);
    }
    ಇಲ್ಲದಿದ್ದರೆ {
      ಹಿಂತಿರುಗಿ ncpu;
    }
  }
}
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *TRANSA, ಅಕ್ಷರ *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t ವಾದಗಳು;

  ಪೂರ್ಣಾಂಕ transa, transb, nrowa, nrowb;
  blasint info;

  ಅಕ್ಷರ transA, transB;
  IFLOAT *ಬಫರ್;
  IFLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ದ್ವಿಗುಣ MNK;
#ಆದರೆ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ಪೂರ್ಣಾಂಕ nodes;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_NAME;

  ವಾದಗಳು.m = *M;
  ವಾದಗಳು.n = *N;
  ವಾದಗಳು.k = *K;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;
  ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

  ವಾದಗಳು.lda = *ldA;
  ವಾದಗಳು.ldb = *ldB;
  ವಾದಗಳು.ldc = *ldC;

  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  ಆದರೆ (transA == 'N') transa = 0;
  ಆದರೆ (transA == 'T') transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ಆದರೆ (transA == 'R') transa = 0;
  ಆದರೆ (transA == 'C') transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ (transA == 'R') transa = 2;
  ಆದರೆ (transA == 'C') transa = 3;
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (transB == 'N') transb = 0;
  ಆದರೆ (transB == 'T') transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ಆದರೆ (transB == 'R') transb = 0;
  ಆದರೆ (transB == 'C') transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ (transB == 'R') transb = 2;
  ಆದರೆ (transB == 'C') transb = 3;
#ಆದರೆ_ಅಂತ್ಯ
  nrowa = ವಾದಗಳು.m;
  ಆದರೆ (transa & 1) nrowa = ವಾದಗಳು.k;
  nrowb = ವಾದಗಳು.k;
  ಆದರೆ (transb & 1) nrowb = ವಾದಗಳು.n;

  info = 0;

  ಆದರೆ (ವಾದಗಳು.ldc < ವಾದಗಳು.m) info = 13;
  ಆದರೆ (ವಾದಗಳು.ldb < nrowb)  info = 10;
  ಆದರೆ (ವಾದಗಳು.lda < nrowa)  info =  8;
  ಆದರೆ (ವಾದಗಳು.k < 0)        info =  5;
  ಆದರೆ (ವಾದಗಳು.n < 0)        info =  4;
  ಆದರೆ (ವಾದಗಳು.m < 0)        info =  3;
  ಆದರೆ (transb < 0)        info =  2;
  ಆದರೆ (transa < 0)        info =  1;

  ಆದರೆ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_TRANSPOSE TransA, ಎಣಿಕೆ CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ಇಲ್ಲದಿದ್ದರೆ
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
#ಆದರೆ_ಅಂತ್ಯ

  blas_arg_t ವಾದಗಳು;
  ಪೂರ್ಣಾಂಕ transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *ಬಫರ್;
  XFLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ದ್ವಿಗುಣ MNK;
#ಆದರೆ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  ಪೂರ್ಣಾಂಕ nodes;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)&alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)&beta;
#ಇಲ್ಲದಿದ್ದರೆ
  ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
  ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;
#ಆದರೆ_ಅಂತ್ಯ

  transa = -1;
  transb = -1;
  info   =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ವಾದಗಳು.m = m;
    ವಾದಗಳು.n = n;
    ವಾದಗಳು.k = k;

    ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
    ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;
    ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

    ವಾದಗಳು.lda = lda;
    ವಾದಗಳು.ldb = ldb;
    ವಾದಗಳು.ldc = ldc;

    ಆದರೆ (TransA == CblasNoTrans)     transa = 0;
    ಆದರೆ (TransA == CblasTrans)       transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (TransA == CblasConjNoTrans) transa = 0;
    ಆದರೆ (TransA == CblasConjTrans)   transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (TransA == CblasConjNoTrans) transa = 2;
    ಆದರೆ (TransA == CblasConjTrans)   transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
    ಆದರೆ (TransB == CblasNoTrans)     transb = 0;
    ಆದರೆ (TransB == CblasTrans)       transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (TransB == CblasConjNoTrans) transb = 0;
    ಆದರೆ (TransB == CblasConjTrans)   transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (TransB == CblasConjNoTrans) transb = 2;
    ಆದರೆ (TransB == CblasConjTrans)   transb = 3;
#ಆದರೆ_ಅಂತ್ಯ

    nrowa = ವಾದಗಳು.m;
    ಆದರೆ (transa & 1) nrowa = ವಾದಗಳು.k;
    nrowb = ವಾದಗಳು.k;
    ಆದರೆ (transb & 1) nrowb = ವಾದಗಳು.n;

    info = -1;

    ಆದರೆ (ವಾದಗಳು.ldc < ವಾದಗಳು.m) info = 13;
    ಆದರೆ (ವಾದಗಳು.ldb < nrowb)  info = 10;
    ಆದರೆ (ವಾದಗಳು.lda < nrowa)  info =  8;
    ಆದರೆ (ವಾದಗಳು.k < 0)        info =  5;
    ಆದರೆ (ವಾದಗಳು.n < 0)        info =  4;
    ಆದರೆ (ವಾದಗಳು.m < 0)        info =  3;
    ಆದರೆ (transb < 0)        info =  2;
    ಆದರೆ (transa < 0)        info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ವಾದಗಳು.m = n;
    ವಾದಗಳು.n = m;
    ವಾದಗಳು.k = k;

    ವಾದಗಳು.a = (ಶೂನ್ಯ *)b;
    ವಾದಗಳು.b = (ಶೂನ್ಯ *)a;
    ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

    ವಾದಗಳು.lda = ldb;
    ವಾದಗಳು.ldb = lda;
    ವಾದಗಳು.ldc = ldc;

    ಆದರೆ (TransB == CblasNoTrans)     transa = 0;
    ಆದರೆ (TransB == CblasTrans)       transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (TransB == CblasConjNoTrans) transa = 0;
    ಆದರೆ (TransB == CblasConjTrans)   transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (TransB == CblasConjNoTrans) transa = 2;
    ಆದರೆ (TransB == CblasConjTrans)   transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
    ಆದರೆ (TransA == CblasNoTrans)     transb = 0;
    ಆದರೆ (TransA == CblasTrans)       transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (TransA == CblasConjNoTrans) transb = 0;
    ಆದರೆ (TransA == CblasConjTrans)   transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (TransA == CblasConjNoTrans) transb = 2;
    ಆದರೆ (TransA == CblasConjTrans)   transb = 3;
#ಆದರೆ_ಅಂತ್ಯ

    nrowa = ವಾದಗಳು.m;
    ಆದರೆ (transa & 1) nrowa = ವಾದಗಳು.k;
    nrowb = ವಾದಗಳು.k;
    ಆದರೆ (transb & 1) nrowb = ವಾದಗಳು.n;

    info = -1;

    ಆದರೆ (ವಾದಗಳು.ldc < ವಾದಗಳು.m) info = 13;
    ಆದರೆ (ವಾದಗಳು.ldb < nrowb)  info = 10;
    ಆದರೆ (ವಾದಗಳು.lda < nrowa)  info =  8;
    ಆದರೆ (ವಾದಗಳು.k < 0)        info =  5;
    ಆದರೆ (ವಾದಗಳು.n < 0)        info =  4;
    ಆದರೆ (ವಾದಗಳು.m < 0)        info =  3;
    ಆದರೆ (transb < 0)        info =  2;
    ಆದರೆ (transa < 0)        info =  1;

  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }


  ಆದರೆ ((ವಾದಗಳು.m == 0) || (ವಾದಗಳು.n == 0)) ಹಿಂತಿರುಗಿ;
#ಆದರೆ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಆದರೆ defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
  ಆದರೆ (support_avx512() )
#ಆದರೆ_ಅಂತ್ಯ
  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ನೇರ_PERFORMANT(m,n,k)) {
        SGEMM_ನೇರ(m, n, k, a, lda, b, ldb, c, ldc);
        ಹಿಂತಿರುಗಿ;
  }
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
ಆದರೆ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಆದರೆ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ಆದರೆ_ಅಂತ್ಯ
)
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ನೇರ_PERFORMANT(m,n,k)) {
        SGEMM_ನೇರ(m, n, k, a, lda, b, ldb, c, ldc);
        ಹಿಂತಿರುಗಿ;
  }
ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ನೇರ_PERFORMANT(m,n,k)) {
        SGEMM_ನೇರ_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ //defined arm64
#ಆದರೆ_ಅಂತ್ಯ //defined complex


#ಆದರೆ_ಅಂತ್ಯ //defined CBLAS

#ಆದರೆ !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಆದರೆ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
ಆದರೆ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಆದರೆ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ಆದರೆ_ಅಂತ್ಯ
)
#ಆದರೆ_ಅಂತ್ಯ //defined dynarch
{
ಅಕ್ಷರ* TA,*TB;
      ಆದರೆ (transa & 1) 
        TA = "T";
	ಇಲ್ಲದಿದ್ದರೆ 
	TA= "N";
      ಆದರೆ (transb & 1) 
        TB = "T";
	ಇಲ್ಲದಿದ್ದರೆ 
	TB= "N";
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transa == 3) 
	TA= "T";
      ಆದರೆ (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)ವಾದಗಳು.alpha;
FLOAT* be=(FLOAT*)ವಾದಗಳು.beta;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ DOUBLE
        SME_SGEMM_KERNEL(TA,TB, ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, al, ವಾದಗಳು.a, ವಾದಗಳು.lda, ವಾದಗಳು.b, ವಾದಗಳು.ldb, be, ವಾದಗಳು.c, ವಾದಗಳು.ldc);
#ಇಲ್ಲದಿದ್ದರೆ
        SME_DGEMM_KERNEL(TA,TB, ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, al, ವಾದಗಳು.a, ವಾದಗಳು.lda, ವಾದಗಳು.b, ವಾದಗಳು.ldb, be, ವಾದಗಳು.c, ವಾದಗಳು.ldc);
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transa == 2) 
	TA= "R";
      ಆದರೆ (transb == 2) 
	TB= "R";
      ಆದರೆ (transa == 3) 
	TA= "C";
      ಆದರೆ (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)ವಾದಗಳು.alpha;
FLOAT* be=(FLOAT*)ವಾದಗಳು.beta;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ DOUBLE
ದಶಮಾಂಶ _Complex c_al={al[0],al[1]};
ದಶಮಾಂಶ _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, al[0],al[1], ವಾದಗಳು.a, ವಾದಗಳು.lda, ವಾದಗಳು.b, ವಾದಗಳು.ldb, be[0],be[1], ವಾದಗಳು.c, ವಾದಗಳು.ldc);
#ಇಲ್ಲದಿದ್ದರೆ
ದ್ವಿಗುಣ _Complex c_al={al[0],al[1]};
ದ್ವಿಗುಣ _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, al[0],al[1], ವಾದಗಳು.a, ವಾದಗಳು.lda, ವಾದಗಳು.b, ವಾದಗಳು.ldb, be[0],be[1], ವಾದಗಳು.c, ವಾದಗಳು.ldc);
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
	ಹಿಂತಿರುಗಿ;
}
#ಆದರೆ_ಅಂತ್ಯ //defined arm64

#ಆದರೆ_ಅಂತ್ಯ //defined b/hfloat16 

#ಆದರೆ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ಆದರೆ defined(DYNAMIC_ARCH)
  ಆದರೆ (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ಹಿಂತಿರುಗಿ;
  }
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  ಆದರೆ (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    ಹಿಂತಿರುಗಿ;
  }
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  ಆದರೆ ((ವಾದಗಳು.m == 0) || (ವಾದಗಳು.n == 0)) ಹಿಂತಿರುಗಿ;

#ಆದರೆ 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, ವಾದಗಳು.lda, ವಾದಗಳು.ldb, ವಾದಗಳು.ldc);
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#ವ್ಯಾಖ್ಯಾನಿಸು BFLOAT16_GEMM_GEMV_FORWARD 1
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು BFLOAT16_GEMM_GEMV_FORWARD 0
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#ವ್ಯಾಖ್ಯಾನಿಸು HFLOAT16_GEMM_GEMV_FORWARD 1
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು HFLOAT16_GEMM_GEMV_FORWARD 0
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#ಆದರೆ defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#ಇಲ್ಲದಿದ್ದರೆ
  bool have_tuned_gemv = true;
#ಆದರೆ_ಅಂತ್ಯ
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  ಸ್ಥಿರ ಸ್ಥಿರಾಂಕ ಅಕ್ಷರ gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  ಆದರೆ (ವಾದಗಳು.k != 0) {
    ಆದರೆ (ವಾದಗಳು.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = ವಾದಗಳು.m;
      blasint n = ವಾದಗಳು.k;
      blasint lda = ವಾದಗಳು.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      ಪೂರ್ಣಾಂಕ trans_idx = transa | ((transb & 2) << 1);
      ಅಕ್ಷರ NT = gemv_trans[trans_idx];
      ಆದರೆ (transa & 1) {
        m = ವಾದಗಳು.k;
        n = ವಾದಗಳು.m;
      }
      ಆದರೆ (transb & 1) {
        inc_x = ವಾದಗಳು.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      ಆದರೆ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, ವಾದಗಳು.alpha, ವಾದಗಳು.a, &lda, ವಾದಗಳು.b, &inc_x, ವಾದಗಳು.beta, ವಾದಗಳು.c, &inc_y);
        ಹಿಂತಿರುಗಿ;
      }
    }
    ಆದರೆ (ವಾದಗಳು.m == 1) {
      blasint inc_x = ವಾದಗಳು.lda;
      blasint inc_y = ವಾದಗಳು.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = ವಾದಗಳು.k;
      blasint n = ವಾದಗಳು.n;
      blasint ldb = ವಾದಗಳು.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      ಪೂರ್ಣಾಂಕ trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      ಅಕ್ಷರ NT = gemv_trans[trans_idx];
      ಆದರೆ (transa & 1) {
        inc_x = 1;
      }
      ಆದರೆ (transb & 1) {
        m = ವಾದಗಳು.n;
        n = ವಾದಗಳು.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      ಆದರೆ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, ವಾದಗಳು.alpha, ವಾದಗಳು.b, &ldb, ವಾದಗಳು.a, &inc_x, ವಾದಗಳು.beta, ವಾದಗಳು.c, &inc_y);
        ಹಿಂತಿರುಗಿ;
      }
    }
  }
#ಆದರೆ_ಅಂತ್ಯ

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

#ಆದರೆ USE_SMALL_ಮಾತೃಕೆ_OPT
ಆದರೆ (*(FLOAT*)ವಾದಗಳು.alpha != 0.0) {
#ಆದರೆ !defined(COMPLEX)
  ಆದರೆ(GEMM_SMALL_ಮಾತೃಕೆ_PERMIT(transa, transb, ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, *(FLOAT *)(ವಾದಗಳು.alpha), *(FLOAT *)(ವಾದಗಳು.beta))){

	  ಆದರೆ(*(FLOAT *)(ವಾದಗಳು.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, ವಾದಗಳು.a, ವಾದಗಳು.lda, *(FLOAT *)(ವಾದಗಳು.alpha), ವಾದಗಳು.b, ವಾದಗಳು.ldb, ವಾದಗಳು.c, ವಾದಗಳು.ldc);
	  }ಇಲ್ಲದಿದ್ದರೆ{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, ವಾದಗಳು.a, ವಾದಗಳು.lda, *(FLOAT *)(ವಾದಗಳು.alpha), ವಾದಗಳು.b, ವಾದಗಳು.ldb, *(FLOAT *)(ವಾದಗಳು.beta), ವಾದಗಳು.c, ವಾದಗಳು.ldc);
	  }
	  ಹಿಂತಿರುಗಿ;
  }
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ(GEMM_SMALL_ಮಾತೃಕೆ_PERMIT(transa, transb, ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, alpha[0], alpha[1], beta[0], beta[1])){
	  ಆದರೆ(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, ವಾದಗಳು.a, ವಾದಗಳು.lda, alpha[0], alpha[1], ವಾದಗಳು.b, ವಾದಗಳು.ldb, ವಾದಗಳು.c, ವಾದಗಳು.ldc);
	  }ಇಲ್ಲದಿದ್ದರೆ{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(ವಾದಗಳು.m, ವಾದಗಳು.n, ವಾದಗಳು.k, ವಾದಗಳು.a, ವಾದಗಳು.lda, alpha[0], alpha[1], ವಾದಗಳು.b, ವಾದಗಳು.ldb, beta[0], beta[1], ವಾದಗಳು.c, ವಾದಗಳು.ldc);
	  }
	  ಹಿಂತಿರುಗಿ;
  }
#ಆದರೆ_ಅಂತ್ಯ
}
#ಆದರೆ_ಅಂತ್ಯ

  ಬಫರ್ = (XFLOAT *)blas_ಸ್ಮರಣೆ_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#ಆದರೆ defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)ಬಫರ್ + (WhereAmI() & 0xf) * GEMM_ಸ್ಥಾನಾಂತರ_A);
#ಇಲ್ಲದಿದ್ದರೆ
  sa = (XFLOAT *)((BLASLONG)ಬಫರ್ +GEMM_ಸ್ಥಾನಾಂತರ_A);
#ಆದರೆ_ಅಂತ್ಯ
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ಗಾತ್ರ_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ಸ್ಥಾನಾಂತರ_B);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ಆದರೆ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  ಬಹುಲಕ |= (transa << BLAS_TRANSA_ಸರಿಸು);
  ಬಹುಲಕ |= (transb << BLAS_TRANSB_ಸರಿಸು);
#ಆದರೆ_ಅಂತ್ಯ

  MNK = (ದ್ವಿಗುಣ) ವಾದಗಳು.m * (ದ್ವಿಗುಣ) ವಾದಗಳು.n * (ದ್ವಿಗುಣ) ವಾದಗಳು.k;
  ವಾದಗಳು.nthreads = get_gemm_optimal_nthreads(MNK);

  ವಾದಗಳು.common = ಖಾಲಿ;

 ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (gemm[(transb << 2) | transa])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP

  } ಇಲ್ಲದಿದ್ದರೆ {

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ USE_SIMPLE_THREADED_LEVEL3

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ NO_AFFINITY
      nodes = get_ಅಂಕೆ_nodes();

      ಆದರೆ ((nodes > 1) && get_ಗಂಟು_equal()) {

	ವಾದಗಳು.nthreads /= nodes;

	gemm_thread_mn(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } ಇಲ್ಲದಿದ್ದರೆ {
#ಆದರೆ_ಅಂತ್ಯ

	(gemm[16 | (transb << 2) | transa])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ಇಲ್ಲದಿದ್ದರೆ

	GEMM_THREAD(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, gemm[(transb << 2) | transa], sa, sb, ವಾದಗಳು.nthreads);

#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ USE_SIMPLE_THREADED_LEVEL3
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ NO_AFFINITY
      }
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  }
#ಆದರೆ_ಅಂತ್ಯ

 blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE, ವಾದಗಳು.m * ವಾದಗಳು.k + ವಾದಗಳು.k * ವಾದಗಳು.n + ವಾದಗಳು.m * ವಾದಗಳು.n, 2 * ವಾದಗಳು.m * ವಾದಗಳು.n * ವಾದಗಳು.k);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
