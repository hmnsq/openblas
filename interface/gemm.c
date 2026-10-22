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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव <stdbool.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_अपरिभाषितम् COMPLEX
#परिभाषय SMP_देहली_MIN 65536.0
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGEMM "
#परिभाषय GEMV BLASFUNC(qgemv)
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGEMM "
#परिभाषय GEMV BLASFUNC(dgemv)
#नोचेद्यदि defined(BFLOAT16)
#यदि_परिभाषितम् BGEMM
#परिभाषय ERROR_NAME "BGEMM "
#परिभाषय GEMV BLASFUNC(bgemv)
#अन्यथा
#परिभाषय ERROR_NAME "SBGEMM "
#परिभाषय GEMV BLASFUNC(sbgemv)
#यदि_अन्तः
#नोचेद्यदि defined(HFLOAT16)
#परिभाषय ERROR_NAME "SHGEMM "
#अन्यथा
#परिभाषय ERROR_NAME "SGEMM "
#परिभाषय GEMV BLASFUNC(sgemv)
#यदि_अन्तः
#अन्यथा
#परिभाषय SMP_देहली_MIN 8192.0
#यदि_अपरिभाषितम् GEMM3M
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMM "
#परिभाषय GEMV BLASFUNC(xgemv)
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMM "
#परिभाषय GEMV BLASFUNC(zgemv)
#अन्यथा
#परिभाषय ERROR_NAME "CGEMM "
#परिभाषय GEMV BLASFUNC(cgemv)
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMM3M "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMM3M "
#अन्यथा
#परिभाषय ERROR_NAME "CGEMM3M "
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि_अपरिभाषितम् GEMM_MULTITHREAD_देहली
#परिभाषय GEMM_MULTITHREAD_देहली 4
#यदि_अन्तः

स्थिर पूर्णाङ्क (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#यदि !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#यदि_अन्तः
#अन्यथा
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#यदि_अन्तः
#यदि_अन्तः
};

#यदि defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#परिभाषय USE_SMALL_MATRIX_OPT 1
#अन्यथा
#परिभाषय USE_SMALL_MATRIX_OPT 0
#यदि_अन्तः

#यदि USE_SMALL_MATRIX_OPT
#यदि_अपरिभाषितम् DYNAMIC_ARCH
#परिभाषय SMALL_KERNEL_ADDR(सारणी, idx) ((शून्य *)(सारणी[idx]))
#अन्यथा
#परिभाषय SMALL_KERNEL_ADDR(सारणी, idx) ((शून्य *)(*(uintptr_t *)((वर्ण *)gotoblas + (size_t)(सारणी[idx]))))
#यदि_अन्तः


#यदि_अपरिभाषितम् COMPLEX
स्थिर size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


स्थिर size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#परिभाषय GEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#परिभाषय GEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

स्थिर size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

स्थिर size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#परिभाषय ZGEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#परिभाषय ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#यदि_अन्तः
#यदि_अन्तः

#यदि defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#परिभाषय XFEATURE_XTILEDATA 18
#परिभाषय ARCH_REQ_XCOMP_PERM 0x1023
स्थिर पूर्णाङ्क openblas_amxtile_permission = 0;
स्थिर पूर्णाङ्क init_amxtile_permission() {
  दीर्घ स्थितिः =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  यदि (स्थितिः != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    प्रत्यागम -1;
  }
  openblas_amxtile_permission = 1;
  प्रत्यागम 0;
}
#यदि_अन्तः

#यदि_परिभाषितम् DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);
#यदि_अन्तः

#यदि_परिभाषितम् SMP
#यदि defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
स्थिर अन्तःस्थ पूर्णाङ्क get_gemm_optimal_nthreads_neoversev1(द्विगुण MNK, पूर्णाङ्क ncpu) {
  प्रत्यागम
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
#यदि_अन्तः

#यदि defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
स्थिर अन्तःस्थ पूर्णाङ्क get_gemm_optimal_nthreads_neoversev2(द्विगुण MNK, पूर्णाङ्क ncpu) {
  प्रत्यागम
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
#यदि_अन्तः

स्थिर अन्तःस्थ पूर्णाङ्क get_gemm_optimal_nthreads(द्विगुण MNK) {
  पूर्णाङ्क ncpu = अङ्कः_cpu_avail(3);
#यदि defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  प्रत्यागम get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#नोचेद्यदि defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  प्रत्यागम get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#नोचेद्यदि defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  यदि (वाक्यं_तुलय(gotoblas_corename(), "neoversev1") == 0) {
    प्रत्यागम get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  यदि (वाक्यं_तुलय(gotoblas_corename(), "neoversev2") == 0) {
    प्रत्यागम get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#यदि_अन्तः
  यदि ( MNK <= (SMP_देहली_MIN  * (द्विगुण) GEMM_MULTITHREAD_देहली) ) {
    प्रत्यागम 1;
  }
  अन्यथा {
    यदि (MNK/ncpu < SMP_देहली_MIN*(द्विगुण)GEMM_MULTITHREAD_देहली) {
      प्रत्यागम MNK/(SMP_देहली_MIN*(द्विगुण)GEMM_MULTITHREAD_देहली);
    }
    अन्यथा {
      प्रत्यागम ncpu;
    }
  }
}
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *TRANSA, वर्ण *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t तर्काः;

  पूर्णाङ्क transa, transb, nrowa, nrowb;
  blasint info;

  वर्ण transA, transB;
  IFLOAT *बफर;
  IFLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
  द्विगुण MNK;
#यदि defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  पूर्णाङ्क nodes;
#यदि_अन्तः

  PRINT_DEBUG_NAME;

  तर्काः.m = *M;
  तर्काः.n = *N;
  तर्काः.k = *K;

  तर्काः.a = (शून्य *)a;
  तर्काः.b = (शून्य *)b;
  तर्काः.c = (शून्य *)c;

  तर्काः.lda = *ldA;
  तर्काः.ldb = *ldB;
  तर्काः.ldc = *ldC;

  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  यदि (transA == 'N') transa = 0;
  यदि (transA == 'T') transa = 1;
#यदि_अपरिभाषितम् COMPLEX
  यदि (transA == 'R') transa = 0;
  यदि (transA == 'C') transa = 1;
#अन्यथा
  यदि (transA == 'R') transa = 2;
  यदि (transA == 'C') transa = 3;
#यदि_अन्तः

  यदि (transB == 'N') transb = 0;
  यदि (transB == 'T') transb = 1;
#यदि_अपरिभाषितम् COMPLEX
  यदि (transB == 'R') transb = 0;
  यदि (transB == 'C') transb = 1;
#अन्यथा
  यदि (transB == 'R') transb = 2;
  यदि (transB == 'C') transb = 3;
#यदि_अन्तः
  nrowa = तर्काः.m;
  यदि (transa & 1) nrowa = तर्काः.k;
  nrowb = तर्काः.k;
  यदि (transb & 1) nrowb = तर्काः.n;

  info = 0;

  यदि (तर्काः.ldc < तर्काः.m) info = 13;
  यदि (तर्काः.ldb < nrowb)  info = 10;
  यदि (तर्काः.lda < nrowa)  info =  8;
  यदि (तर्काः.k < 0)        info =  5;
  यदि (तर्काः.n < 0)        info =  4;
  यदि (तर्काः.m < 0)        info =  3;
  यदि (transb < 0)        info =  2;
  यदि (transa < 0)        info =  1;

  यदि (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb,
	   शून्य *vbeta,
	   शून्य *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#यदि_अन्तः

  blas_arg_t तर्काः;
  पूर्णाङ्क transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *बफर;
  XFLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
  द्विगुण MNK;
#यदि defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  पूर्णाङ्क nodes;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

#यदि_अपरिभाषितम् COMPLEX
  तर्काः.alpha = (शून्य *)&alpha;
  तर्काः.beta  = (शून्य *)&beta;
#अन्यथा
  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;
#यदि_अन्तः

  transa = -1;
  transb = -1;
  info   =  0;

  यदि (क्रमः == CblasColMajor) {
    तर्काः.m = m;
    तर्काः.n = n;
    तर्काः.k = k;

    तर्काः.a = (शून्य *)a;
    तर्काः.b = (शून्य *)b;
    तर्काः.c = (शून्य *)c;

    तर्काः.lda = lda;
    तर्काः.ldb = ldb;
    तर्काः.ldc = ldc;

    यदि (TransA == CblasNoTrans)     transa = 0;
    यदि (TransA == CblasTrans)       transa = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (TransA == CblasConjNoTrans) transa = 0;
    यदि (TransA == CblasConjTrans)   transa = 1;
#अन्यथा
    यदि (TransA == CblasConjNoTrans) transa = 2;
    यदि (TransA == CblasConjTrans)   transa = 3;
#यदि_अन्तः
    यदि (TransB == CblasNoTrans)     transb = 0;
    यदि (TransB == CblasTrans)       transb = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (TransB == CblasConjNoTrans) transb = 0;
    यदि (TransB == CblasConjTrans)   transb = 1;
#अन्यथा
    यदि (TransB == CblasConjNoTrans) transb = 2;
    यदि (TransB == CblasConjTrans)   transb = 3;
#यदि_अन्तः

    nrowa = तर्काः.m;
    यदि (transa & 1) nrowa = तर्काः.k;
    nrowb = तर्काः.k;
    यदि (transb & 1) nrowb = तर्काः.n;

    info = -1;

    यदि (तर्काः.ldc < तर्काः.m) info = 13;
    यदि (तर्काः.ldb < nrowb)  info = 10;
    यदि (तर्काः.lda < nrowa)  info =  8;
    यदि (तर्काः.k < 0)        info =  5;
    यदि (तर्काः.n < 0)        info =  4;
    यदि (तर्काः.m < 0)        info =  3;
    यदि (transb < 0)        info =  2;
    यदि (transa < 0)        info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {
    तर्काः.m = n;
    तर्काः.n = m;
    तर्काः.k = k;

    तर्काः.a = (शून्य *)b;
    तर्काः.b = (शून्य *)a;
    तर्काः.c = (शून्य *)c;

    तर्काः.lda = ldb;
    तर्काः.ldb = lda;
    तर्काः.ldc = ldc;

    यदि (TransB == CblasNoTrans)     transa = 0;
    यदि (TransB == CblasTrans)       transa = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (TransB == CblasConjNoTrans) transa = 0;
    यदि (TransB == CblasConjTrans)   transa = 1;
#अन्यथा
    यदि (TransB == CblasConjNoTrans) transa = 2;
    यदि (TransB == CblasConjTrans)   transa = 3;
#यदि_अन्तः
    यदि (TransA == CblasNoTrans)     transb = 0;
    यदि (TransA == CblasTrans)       transb = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (TransA == CblasConjNoTrans) transb = 0;
    यदि (TransA == CblasConjTrans)   transb = 1;
#अन्यथा
    यदि (TransA == CblasConjNoTrans) transb = 2;
    यदि (TransA == CblasConjTrans)   transb = 3;
#यदि_अन्तः

    nrowa = तर्काः.m;
    यदि (transa & 1) nrowa = तर्काः.k;
    nrowb = तर्काः.k;
    यदि (transb & 1) nrowb = तर्काः.n;

    info = -1;

    यदि (तर्काः.ldc < तर्काः.m) info = 13;
    यदि (तर्काः.ldb < nrowb)  info = 10;
    यदि (तर्काः.lda < nrowa)  info =  8;
    यदि (तर्काः.k < 0)        info =  5;
    यदि (तर्काः.n < 0)        info =  4;
    यदि (तर्काः.m < 0)        info =  3;
    यदि (transb < 0)        info =  2;
    यदि (transa < 0)        info =  1;

  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }


  यदि ((तर्काः.m == 0) || (तर्काः.n == 0)) प्रत्यागम;
#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
  यदि (support_avx512() )
#यदि_अन्तः
  यदि (क्रमः == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_प्रत्यक्षम्_PERFORMANT(m,n,k)) {
        SGEMM_प्रत्यक्षम्(m, n, k, a, lda, b, ldb, c, ldc);
        प्रत्यागम;
  }
#यदि_अन्तः
#यदि defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (वाक्यं_तुलय(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || वाक्यं_तुलय(gotoblas_corename(), "vortexm4") == 0
#यदि_अन्तः
)
#यदि_अन्तः

  यदि (क्रमः == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_प्रत्यक्षम्_PERFORMANT(m,n,k)) {
        SGEMM_प्रत्यक्षम्(m, n, k, a, lda, b, ldb, c, ldc);
        प्रत्यागम;
  }
अन्यथा यदि (क्रमः == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_प्रत्यक्षम्_PERFORMANT(m,n,k)) {
        SGEMM_प्रत्यक्षम्_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
  }

#यदि_अन्तः //defined arm64
#यदि_अन्तः //defined complex


#यदि_अन्तः //defined CBLAS

#यदि !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (वाक्यं_तुलय(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || वाक्यं_तुलय(gotoblas_corename(), "vortexm4") == 0
#यदि_अन्तः
)
#यदि_अन्तः //defined dynarch
{
वर्ण* TA,*TB;
      यदि (transa & 1) 
        TA = "T";
	अन्यथा 
	TA= "N";
      यदि (transb & 1) 
        TB = "T";
	अन्यथा 
	TB= "N";
#यदि_अपरिभाषितम् COMPLEX
      यदि (transa == 3) 
	TA= "T";
      यदि (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)तर्काः.alpha;
FLOAT* be=(FLOAT*)तर्काः.beta;
#यदि_अपरिभाषितम् DOUBLE
        SME_SGEMM_KERNEL(TA,TB, तर्काः.m, तर्काः.n, तर्काः.k, al, तर्काः.a, तर्काः.lda, तर्काः.b, तर्काः.ldb, be, तर्काः.c, तर्काः.ldc);
#अन्यथा
        SME_DGEMM_KERNEL(TA,TB, तर्काः.m, तर्काः.n, तर्काः.k, al, तर्काः.a, तर्काः.lda, तर्काः.b, तर्काः.ldb, be, तर्काः.c, तर्काः.ldc);
#यदि_अन्तः
#अन्यथा
      यदि (transa == 2) 
	TA= "R";
      यदि (transb == 2) 
	TB= "R";
      यदि (transa == 3) 
	TA= "C";
      यदि (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)तर्काः.alpha;
FLOAT* be=(FLOAT*)तर्काः.beta;
#यदि_अपरिभाषितम् DOUBLE
भिन्न _Complex c_al={al[0],al[1]};
भिन्न _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, तर्काः.m, तर्काः.n, तर्काः.k, al[0],al[1], तर्काः.a, तर्काः.lda, तर्काः.b, तर्काः.ldb, be[0],be[1], तर्काः.c, तर्काः.ldc);
#अन्यथा
द्विगुण _Complex c_al={al[0],al[1]};
द्विगुण _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, तर्काः.m, तर्काः.n, तर्काः.k, al[0],al[1], तर्काः.a, तर्काः.lda, तर्काः.b, तर्काः.ldb, be[0],be[1], तर्काः.c, तर्काः.ldc);
#यदि_अन्तः
#यदि_अन्तः
	प्रत्यागम;
}
#यदि_अन्तः //defined arm64

#यदि_अन्तः //defined b/hfloat16 

#यदि defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#यदि defined(DYNAMIC_ARCH)
  यदि (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    प्रत्यागम;
  }
#यदि_अन्तः
#यदि !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  यदि (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    प्रत्यागम;
  }
#यदि_अन्तः
#यदि_अन्तः  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  यदि ((तर्काः.m == 0) || (तर्काः.n == 0)) प्रत्यागम;

#यदि 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 तर्काः.m, तर्काः.n, तर्काः.k, तर्काः.lda, तर्काः.ldb, तर्काः.ldc);
#यदि_अन्तः

#यदि (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#परिभाषय BFLOAT16_GEMM_GEMV_FORWARD 1
#अन्यथा
#परिभाषय BFLOAT16_GEMM_GEMV_FORWARD 0
#यदि_अन्तः
#यदि (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#परिभाषय HFLOAT16_GEMM_GEMV_FORWARD 1
#अन्यथा
#परिभाषय HFLOAT16_GEMM_GEMV_FORWARD 0
#यदि_अन्तः

#यदि defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#यदि defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#अन्यथा
  bool have_tuned_gemv = true;
#यदि_अन्तः
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  स्थिर ध्रुव वर्ण gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  यदि (तर्काः.k != 0) {
    यदि (तर्काः.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = तर्काः.m;
      blasint n = तर्काः.k;
      blasint lda = तर्काः.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      पूर्णाङ्क trans_idx = transa | ((transb & 2) << 1);
      वर्ण NT = gemv_trans[trans_idx];
      यदि (transa & 1) {
        m = तर्काः.k;
        n = तर्काः.m;
      }
      यदि (transb & 1) {
        inc_x = तर्काः.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      यदि (is_efficient_gemv) {
        GEMV(&NT, &m, &n, तर्काः.alpha, तर्काः.a, &lda, तर्काः.b, &inc_x, तर्काः.beta, तर्काः.c, &inc_y);
        प्रत्यागम;
      }
    }
    यदि (तर्काः.m == 1) {
      blasint inc_x = तर्काः.lda;
      blasint inc_y = तर्काः.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = तर्काः.k;
      blasint n = तर्काः.n;
      blasint ldb = तर्काः.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      पूर्णाङ्क trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      वर्ण NT = gemv_trans[trans_idx];
      यदि (transa & 1) {
        inc_x = 1;
      }
      यदि (transb & 1) {
        m = तर्काः.n;
        n = तर्काः.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      यदि (is_efficient_gemv) {
        GEMV(&NT, &m, &n, तर्काः.alpha, तर्काः.b, &ldb, तर्काः.a, &inc_x, तर्काः.beta, तर्काः.c, &inc_y);
        प्रत्यागम;
      }
    }
  }
#यदि_अन्तः

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#यदि USE_SMALL_MATRIX_OPT
यदि (*(FLOAT*)तर्काः.alpha != 0.0) {
#यदि !defined(COMPLEX)
  यदि(GEMM_SMALL_MATRIX_PERMIT(transa, transb, तर्काः.m, तर्काः.n, तर्काः.k, *(FLOAT *)(तर्काः.alpha), *(FLOAT *)(तर्काः.beta))){

	  यदि(*(FLOAT *)(तर्काः.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(तर्काः.m, तर्काः.n, तर्काः.k, तर्काः.a, तर्काः.lda, *(FLOAT *)(तर्काः.alpha), तर्काः.b, तर्काः.ldb, तर्काः.c, तर्काः.ldc);
	  }अन्यथा{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(तर्काः.m, तर्काः.n, तर्काः.k, तर्काः.a, तर्काः.lda, *(FLOAT *)(तर्काः.alpha), तर्काः.b, तर्काः.ldb, *(FLOAT *)(तर्काः.beta), तर्काः.c, तर्काः.ldc);
	  }
	  प्रत्यागम;
  }
#अन्यथा
  यदि(GEMM_SMALL_MATRIX_PERMIT(transa, transb, तर्काः.m, तर्काः.n, तर्काः.k, alpha[0], alpha[1], beta[0], beta[1])){
	  यदि(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(तर्काः.m, तर्काः.n, तर्काः.k, तर्काः.a, तर्काः.lda, alpha[0], alpha[1], तर्काः.b, तर्काः.ldb, तर्काः.c, तर्काः.ldc);
	  }अन्यथा{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(तर्काः.m, तर्काः.n, तर्काः.k, तर्काः.a, तर्काः.lda, alpha[0], alpha[1], तर्काः.b, तर्काः.ldb, beta[0], beta[1], तर्काः.c, तर्काः.ldc);
	  }
	  प्रत्यागम;
  }
#यदि_अन्तः
}
#यदि_अन्तः

  बफर = (XFLOAT *)blas_स्मृतिः_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#यदि defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)बफर + (WhereAmI() & 0xf) * GEMM_विस्थापनम्_A);
#अन्यथा
  sa = (XFLOAT *)((BLASLONG)बफर +GEMM_विस्थापनम्_A);
#यदि_अन्तः
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);

#यदि_परिभाषितम् SMP
#यदि defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  बहुलकम् |= (transa << BLAS_TRANSA_सरणम्);
  बहुलकम् |= (transb << BLAS_TRANSB_सरणम्);
#यदि_अन्तः

  MNK = (द्विगुण) तर्काः.m * (द्विगुण) तर्काः.n * (द्विगुण) तर्काः.k;
  तर्काः.nthreads = get_gemm_optimal_nthreads(MNK);

  तर्काः.common = रिक्त;

 यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

    (gemm[(transb << 2) | transa])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP

  } अन्यथा {

#यदि_अपरिभाषितम् USE_SIMPLE_THREADED_LEVEL3

#यदि_अपरिभाषितम् NO_AFFINITY
      nodes = get_अङ्कः_nodes();

      यदि ((nodes > 1) && get_ग्रन्थिः_equal()) {

	तर्काः.nthreads /= nodes;

	gemm_thread_mn(बहुलकम्, &तर्काः, रिक्त, रिक्त, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } अन्यथा {
#यदि_अन्तः

	(gemm[16 | (transb << 2) | transa])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#अन्यथा

	GEMM_THREAD(बहुलकम्, &तर्काः, रिक्त, रिक्त, gemm[(transb << 2) | transa], sa, sb, तर्काः.nthreads);

#यदि_अन्तः

#यदि_अपरिभाषितम् USE_SIMPLE_THREADED_LEVEL3
#यदि_अपरिभाषितम् NO_AFFINITY
      }
#यदि_अन्तः
#यदि_अन्तः

#यदि_अन्तः

#यदि_परिभाषितम् SMP
  }
#यदि_अन्तः

 blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, तर्काः.m * तर्काः.k + तर्काः.k * तर्काः.n + तर्काः.m * तर्काः.n, 2 * तर्काः.m * तर्काः.n * तर्काः.k);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
