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
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#define GEMV BLASFUNC(qgemv)
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#define GEMV BLASFUNC(dgemv)
#किंवाजर defined(BFLOAT16)
#ifdef BGEMM
#define ERROR_NAME "BGEMM "
#define GEMV BLASFUNC(bgemv)
#अन्यथा
#define ERROR_NAME "SBGEMM "
#define GEMV BLASFUNC(sbgemv)
#endif
#किंवाजर defined(HFLOAT16)
#define ERROR_NAME "SHGEMM "
#अन्यथा
#define ERROR_NAME "SGEMM "
#define GEMV BLASFUNC(sgemv)
#endif
#अन्यथा
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#define GEMV BLASFUNC(xgemv)
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#define GEMV BLASFUNC(zgemv)
#अन्यथा
#define ERROR_NAME "CGEMM "
#define GEMV BLASFUNC(cgemv)
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#अन्यथा
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

स्थिर पूर्णांक (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#जर !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#जर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#endif
#अन्यथा
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#जर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#endif
#endif
};

#जर defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#define USE_SMALL_MATRIX_OPT 1
#अन्यथा
#define USE_SMALL_MATRIX_OPT 0
#endif

#जर USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(table[idx]))
#अन्यथा
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(*(uintptr_t *)((अक्षर *)gotoblas + (माप_प्रकार)(table[idx]))))
#endif


#ifndef COMPLEX
स्थिर माप_प्रकार gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


स्थिर माप_प्रकार gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

स्थिर माप_प्रकार zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

स्थिर माप_प्रकार zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

#जर defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#define XFEATURE_XTILEDATA 18
#define ARCH_REQ_XCOMP_PERM 0x1023
स्थिर पूर्णांक openblas_amxtile_permission = 0;
स्थिर पूर्णांक init_amxtile_permission() {
  दीर्घ status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  जर (status != 0) {
    ख_लिखो(मानक_त्रुटि, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    परत -1;
  }
  openblas_amxtile_permission = 1;
  परत 0;
}
#endif

#ifdef DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(शून्य);
#endif

#ifdef SMP
#जर defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
स्थिर ओळीत पूर्णांक get_gemm_optimal_nthreads_neoversev1(द्विगुण MNK, पूर्णांक ncpu) {
  परत
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

#जर defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
स्थिर ओळीत पूर्णांक get_gemm_optimal_nthreads_neoversev2(द्विगुण MNK, पूर्णांक ncpu) {
  परत
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

स्थिर ओळीत पूर्णांक get_gemm_optimal_nthreads(द्विगुण MNK) {
  पूर्णांक ncpu = अंक_cpu_avail(3);
#जर defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  परत get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#किंवाजर defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  परत get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#किंवाजर defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  जर (म_भेद(gotoblas_corename(), "neoversev1") == 0) {
    परत get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  जर (म_भेद(gotoblas_corename(), "neoversev2") == 0) {
    परत get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#endif
  जर ( MNK <= (SMP_THRESHOLD_MIN  * (द्विगुण) GEMM_MULTITHREAD_THRESHOLD) ) {
    परत 1;
  }
  अन्यथा {
    जर (MNK/ncpu < SMP_THRESHOLD_MIN*(द्विगुण)GEMM_MULTITHREAD_THRESHOLD) {
      परत MNK/(SMP_THRESHOLD_MIN*(द्विगुण)GEMM_MULTITHREAD_THRESHOLD);
    }
    अन्यथा {
      परत ncpu;
    }
  }
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *TRANSA, अक्षर *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t args;

  पूर्णांक transa, transb, nrowa, nrowb;
  blasint info;

  अक्षर transA, transB;
  IFLOAT *buffer;
  IFLOAT *sa, *sb;

#ifdef SMP
  द्विगुण MNK;
#जर defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#जर defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  पूर्णांक nodes;
#endif

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (शून्य *)a;
  args.b = (शून्य *)b;
  args.c = (शून्य *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  जर (transA == 'N') transa = 0;
  जर (transA == 'T') transa = 1;
#ifndef COMPLEX
  जर (transA == 'R') transa = 0;
  जर (transA == 'C') transa = 1;
#अन्यथा
  जर (transA == 'R') transa = 2;
  जर (transA == 'C') transa = 3;
#endif

  जर (transB == 'N') transb = 0;
  जर (transB == 'T') transb = 1;
#ifndef COMPLEX
  जर (transB == 'R') transb = 0;
  जर (transB == 'C') transb = 1;
#अन्यथा
  जर (transB == 'R') transb = 2;
  जर (transB == 'C') transb = 3;
#endif
  nrowa = args.m;
  जर (transa & 1) nrowa = args.k;
  nrowb = args.k;
  जर (transb & 1) nrowb = args.n;

  info = 0;

  जर (args.ldc < args.m) info = 13;
  जर (args.ldb < nrowb)  info = 10;
  जर (args.lda < nrowa)  info =  8;
  जर (args.k < 0)        info =  5;
  जर (args.n < 0)        info =  4;
  जर (args.m < 0)        info =  3;
  जर (transb < 0)        info =  2;
  जर (transa < 0)        info =  1;

  जर (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
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
#endif

  blas_arg_t args;
  पूर्णांक transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  द्विगुण MNK;
#जर defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

#जर defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  पूर्णांक nodes;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (शून्य *)&alpha;
  args.beta  = (शून्य *)&beta;
#अन्यथा
  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;
#endif

  transa = -1;
  transb = -1;
  info   =  0;

  जर (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (शून्य *)a;
    args.b = (शून्य *)b;
    args.c = (शून्य *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    जर (TransA == CblasNoTrans)     transa = 0;
    जर (TransA == CblasTrans)       transa = 1;
#ifndef COMPLEX
    जर (TransA == CblasConjNoTrans) transa = 0;
    जर (TransA == CblasConjTrans)   transa = 1;
#अन्यथा
    जर (TransA == CblasConjNoTrans) transa = 2;
    जर (TransA == CblasConjTrans)   transa = 3;
#endif
    जर (TransB == CblasNoTrans)     transb = 0;
    जर (TransB == CblasTrans)       transb = 1;
#ifndef COMPLEX
    जर (TransB == CblasConjNoTrans) transb = 0;
    जर (TransB == CblasConjTrans)   transb = 1;
#अन्यथा
    जर (TransB == CblasConjNoTrans) transb = 2;
    जर (TransB == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    जर (transa & 1) nrowa = args.k;
    nrowb = args.k;
    जर (transb & 1) nrowb = args.n;

    info = -1;

    जर (args.ldc < args.m) info = 13;
    जर (args.ldb < nrowb)  info = 10;
    जर (args.lda < nrowa)  info =  8;
    जर (args.k < 0)        info =  5;
    जर (args.n < 0)        info =  4;
    जर (args.m < 0)        info =  3;
    जर (transb < 0)        info =  2;
    जर (transa < 0)        info =  1;
  }

  जर (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (शून्य *)b;
    args.b = (शून्य *)a;
    args.c = (शून्य *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    जर (TransB == CblasNoTrans)     transa = 0;
    जर (TransB == CblasTrans)       transa = 1;
#ifndef COMPLEX
    जर (TransB == CblasConjNoTrans) transa = 0;
    जर (TransB == CblasConjTrans)   transa = 1;
#अन्यथा
    जर (TransB == CblasConjNoTrans) transa = 2;
    जर (TransB == CblasConjTrans)   transa = 3;
#endif
    जर (TransA == CblasNoTrans)     transb = 0;
    जर (TransA == CblasTrans)       transb = 1;
#ifndef COMPLEX
    जर (TransA == CblasConjNoTrans) transb = 0;
    जर (TransA == CblasConjTrans)   transb = 1;
#अन्यथा
    जर (TransA == CblasConjNoTrans) transb = 2;
    जर (TransA == CblasConjTrans)   transb = 3;
#endif

    nrowa = args.m;
    जर (transa & 1) nrowa = args.k;
    nrowb = args.k;
    जर (transb & 1) nrowb = args.n;

    info = -1;

    जर (args.ldc < args.m) info = 13;
    जर (args.ldb < nrowb)  info = 10;
    जर (args.lda < nrowa)  info =  8;
    जर (args.k < 0)        info =  5;
    जर (args.n < 0)        info =  4;
    जर (args.m < 0)        info =  3;
    जर (transb < 0)        info =  2;
    जर (transa < 0)        info =  1;

  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }


  जर ((args.m == 0) || (args.n == 0)) परत;
#जर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#जर defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जर defined(DYNAMIC_ARCH)
  जर (support_avx512() )
#endif
  जर (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        परत;
  }
#endif
#जर defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जर defined(DYNAMIC_ARCH)
जर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#जर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif

  जर (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        परत;
  }
अन्यथा जर (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        परत;
  }

#endif //defined arm64
#endif //defined complex


#endif //defined CBLAS

#जर !defined(BFLOAT16)  && !defined(HFLOAT16)
#जर defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जर defined(DYNAMIC_ARCH)
जर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#जर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif //defined dynarch
{
अक्षर* TA,*TB;
      जर (transa & 1) 
        TA = "T";
	अन्यथा 
	TA= "N";
      जर (transb & 1) 
        TB = "T";
	अन्यथा 
	TB= "N";
#ifndef COMPLEX
      जर (transa == 3) 
	TA= "T";
      जर (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#अन्यथा
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#endif
#अन्यथा
      जर (transa == 2) 
	TA= "R";
      जर (transb == 2) 
	TB= "R";
      जर (transa == 3) 
	TA= "C";
      जर (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#ifndef DOUBLE
दशांश _Complex c_al={al[0],al[1]};
दशांश _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#अन्यथा
द्विगुण _Complex c_al={al[0],al[1]};
द्विगुण _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#endif
#endif
	परत;
}
#endif //defined arm64

#endif //defined b/hfloat16 

#जर defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#जर defined(DYNAMIC_ARCH)
  जर (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    परत;
  }
#endif
#जर !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  जर (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    परत;
  }
#endif
#endif  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  जर ((args.m == 0) || (args.n == 0)) परत;

#जर 0
  ख_लिखो(मानक_त्रुटि, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#endif

#जर (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#define BFLOAT16_GEMM_GEMV_FORWARD 1
#अन्यथा
#define BFLOAT16_GEMM_GEMV_FORWARD 0
#endif
#जर (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#define HFLOAT16_GEMM_GEMV_FORWARD 1
#अन्यथा
#define HFLOAT16_GEMM_GEMV_FORWARD 0
#endif

#जर defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#जर defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#अन्यथा
  bool have_tuned_gemv = true;
#endif
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  स्थिर स्थिरांक अक्षर gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  जर (args.k != 0) {
    जर (args.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.m;
      blasint n = args.k;
      blasint lda = args.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      पूर्णांक trans_idx = transa | ((transb & 2) << 1);
      अक्षर NT = gemv_trans[trans_idx];
      जर (transa & 1) {
        m = args.k;
        n = args.m;
      }
      जर (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      जर (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        परत;
      }
    }
    जर (args.m == 1) {
      blasint inc_x = args.lda;
      blasint inc_y = args.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = args.k;
      blasint n = args.n;
      blasint ldb = args.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      पूर्णांक trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      अक्षर NT = gemv_trans[trans_idx];
      जर (transa & 1) {
        inc_x = 1;
      }
      जर (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      जर (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        परत;
      }
    }
  }
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#जर USE_SMALL_MATRIX_OPT
जर (*(FLOAT*)args.alpha != 0.0) {
#जर !defined(COMPLEX)
  जर(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  जर(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }अन्यथा{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  परत;
  }
#अन्यथा
  जर(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  जर(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }अन्यथा{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  परत;
  }
#endif
}
#endif

  buffer = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#जर defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)buffer + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#अन्यथा
  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
#endif
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#जर defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#endif

  MNK = (द्विगुण) args.m * (द्विगुण) args.n * (द्विगुण) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = रिक्त;

 जर (args.nthreads == 1) {
#endif

    (gemm[(transb << 2) | transa])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP

  } अन्यथा {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_अंक_nodes();

      जर ((nodes > 1) && get_गाठ_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, रिक्त, रिक्त, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } अन्यथा {
#endif

	(gemm[16 | (transb << 2) | transa])(&args, रिक्त, रिक्त, sa, sb, 0);

#अन्यथा

	GEMM_THREAD(mode, &args, रिक्त, रिक्त, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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

  परत;
}
