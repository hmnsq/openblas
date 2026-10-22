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

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश <stdbool.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_अघोषित COMPLEX
#घोषणा SMP_THRESHOLD_MIN 65536.0
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGEMM "
#घोषणा GEMV BLASFUNC(qgemv)
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGEMM "
#घोषणा GEMV BLASFUNC(dgemv)
#या_अगर defined(BFLOAT16)
#अगर_घोषित BGEMM
#घोषणा ERROR_NAME "BGEMM "
#घोषणा GEMV BLASFUNC(bgemv)
#अन्यथा
#घोषणा ERROR_NAME "SBGEMM "
#घोषणा GEMV BLASFUNC(sbgemv)
#पूर्ण_अगर
#या_अगर defined(HFLOAT16)
#घोषणा ERROR_NAME "SHGEMM "
#अन्यथा
#घोषणा ERROR_NAME "SGEMM "
#घोषणा GEMV BLASFUNC(sgemv)
#पूर्ण_अगर
#अन्यथा
#घोषणा SMP_THRESHOLD_MIN 8192.0
#अगर_अघोषित GEMM3M
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XGEMM "
#घोषणा GEMV BLASFUNC(xgemv)
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMM "
#घोषणा GEMV BLASFUNC(zgemv)
#अन्यथा
#घोषणा ERROR_NAME "CGEMM "
#घोषणा GEMV BLASFUNC(cgemv)
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XGEMM3M "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMM3M "
#अन्यथा
#घोषणा ERROR_NAME "CGEMM3M "
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर_अघोषित GEMM_MULTITHREAD_THRESHOLD
#घोषणा GEMM_MULTITHREAD_THRESHOLD 4
#पूर्ण_अगर

अटल पूर्णांक (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#अगर !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#पूर्ण_अगर
#अन्यथा
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#पूर्ण_अगर
#पूर्ण_अगर
};

#अगर defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#घोषणा USE_SMALL_MATRIX_OPT 1
#अन्यथा
#घोषणा USE_SMALL_MATRIX_OPT 0
#पूर्ण_अगर

#अगर USE_SMALL_MATRIX_OPT
#अगर_अघोषित DYNAMIC_ARCH
#घोषणा SMALL_KERNEL_ADDR(table, idx) ((व्योम *)(table[idx]))
#अन्यथा
#घोषणा SMALL_KERNEL_ADDR(table, idx) ((व्योम *)(*(uintptr_t *)((अक्षर *)gotoblas + (माप_प्रकार)(table[idx]))))
#पूर्ण_अगर


#अगर_अघोषित COMPLEX
अटल माप_प्रकार gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


अटल माप_प्रकार gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#घोषणा GEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#घोषणा GEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

अटल माप_प्रकार zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

अटल माप_प्रकार zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#घोषणा ZGEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#घोषणा ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#पूर्ण_अगर
#पूर्ण_अगर

#अगर defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#घोषणा XFEATURE_XTILEDATA 18
#घोषणा ARCH_REQ_XCOMP_PERM 0x1023
अटल पूर्णांक openblas_amxtile_permission = 0;
अटल पूर्णांक init_amxtile_permission() {
  दीर्घ status =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  अगर (status != 0) {
    ख_लिखो(मानक_त्रुटि, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    वापस -1;
  }
  openblas_amxtile_permission = 1;
  वापस 0;
}
#पूर्ण_अगर

#अगर_घोषित DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(व्योम);
#पूर्ण_अगर

#अगर_घोषित SMP
#अगर defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
अटल अंतरभूत पूर्णांक get_gemm_optimal_nthreads_neoversev1(द्विगुन MNK, पूर्णांक ncpu) {
  वापस
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
#पूर्ण_अगर

#अगर defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
अटल अंतरभूत पूर्णांक get_gemm_optimal_nthreads_neoversev2(द्विगुन MNK, पूर्णांक ncpu) {
  वापस
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
#पूर्ण_अगर

अटल अंतरभूत पूर्णांक get_gemm_optimal_nthreads(द्विगुन MNK) {
  पूर्णांक ncpu = अंक_cpu_avail(3);
#अगर defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  वापस get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#या_अगर defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  वापस get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#या_अगर defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  अगर (म_भेद(gotoblas_corename(), "neoversev1") == 0) {
    वापस get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  अगर (म_भेद(gotoblas_corename(), "neoversev2") == 0) {
    वापस get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#पूर्ण_अगर
  अगर ( MNK <= (SMP_THRESHOLD_MIN  * (द्विगुन) GEMM_MULTITHREAD_THRESHOLD) ) {
    वापस 1;
  }
  अन्यथा {
    अगर (MNK/ncpu < SMP_THRESHOLD_MIN*(द्विगुन)GEMM_MULTITHREAD_THRESHOLD) {
      वापस MNK/(SMP_THRESHOLD_MIN*(द्विगुन)GEMM_MULTITHREAD_THRESHOLD);
    }
    अन्यथा {
      वापस ncpu;
    }
  }
}
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *TRANSA, अक्षर *TRANSB,
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
  IFLOAT *बफर;
  IFLOAT *sa, *sb;

#अगर_घोषित SMP
  द्विगुन MNK;
#अगर defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  पूर्णांक nodes;
#पूर्ण_अगर

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;
  args.k = *K;

  args.a = (व्योम *)a;
  args.b = (व्योम *)b;
  args.c = (व्योम *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (व्योम *)alpha;
  args.beta  = (व्योम *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  अगर (transA == 'N') transa = 0;
  अगर (transA == 'T') transa = 1;
#अगर_अघोषित COMPLEX
  अगर (transA == 'R') transa = 0;
  अगर (transA == 'C') transa = 1;
#अन्यथा
  अगर (transA == 'R') transa = 2;
  अगर (transA == 'C') transa = 3;
#पूर्ण_अगर

  अगर (transB == 'N') transb = 0;
  अगर (transB == 'T') transb = 1;
#अगर_अघोषित COMPLEX
  अगर (transB == 'R') transb = 0;
  अगर (transB == 'C') transb = 1;
#अन्यथा
  अगर (transB == 'R') transb = 2;
  अगर (transB == 'C') transb = 3;
#पूर्ण_अगर
  nrowa = args.m;
  अगर (transa & 1) nrowa = args.k;
  nrowb = args.k;
  अगर (transb & 1) nrowb = args.n;

  info = 0;

  अगर (args.ldc < args.m) info = 13;
  अगर (args.ldb < nrowb)  info = 10;
  अगर (args.lda < nrowa)  info =  8;
  अगर (args.k < 0)        info =  5;
  अगर (args.n < 0)        info =  4;
  अगर (args.m < 0)        info =  3;
  अगर (transb < 0)        info =  2;
  अगर (transa < 0)        info =  1;

  अगर (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_TRANSPOSE TransA, क्रमागत CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#अगर_अघोषित COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   व्योम *valpha,
	   व्योम *va, blasint lda,
	   व्योम *vb, blasint ldb,
	   व्योम *vbeta,
	   व्योम *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#पूर्ण_अगर

  blas_arg_t args;
  पूर्णांक transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *बफर;
  XFLOAT *sa, *sb;

#अगर_घोषित SMP
  द्विगुन MNK;
#अगर defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  पूर्णांक nodes;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

#अगर_अघोषित COMPLEX
  args.alpha = (व्योम *)&alpha;
  args.beta  = (व्योम *)&beta;
#अन्यथा
  args.alpha = (व्योम *)alpha;
  args.beta  = (व्योम *)beta;
#पूर्ण_अगर

  transa = -1;
  transb = -1;
  info   =  0;

  अगर (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (व्योम *)a;
    args.b = (व्योम *)b;
    args.c = (व्योम *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    अगर (TransA == CblasNoTrans)     transa = 0;
    अगर (TransA == CblasTrans)       transa = 1;
#अगर_अघोषित COMPLEX
    अगर (TransA == CblasConjNoTrans) transa = 0;
    अगर (TransA == CblasConjTrans)   transa = 1;
#अन्यथा
    अगर (TransA == CblasConjNoTrans) transa = 2;
    अगर (TransA == CblasConjTrans)   transa = 3;
#पूर्ण_अगर
    अगर (TransB == CblasNoTrans)     transb = 0;
    अगर (TransB == CblasTrans)       transb = 1;
#अगर_अघोषित COMPLEX
    अगर (TransB == CblasConjNoTrans) transb = 0;
    अगर (TransB == CblasConjTrans)   transb = 1;
#अन्यथा
    अगर (TransB == CblasConjNoTrans) transb = 2;
    अगर (TransB == CblasConjTrans)   transb = 3;
#पूर्ण_अगर

    nrowa = args.m;
    अगर (transa & 1) nrowa = args.k;
    nrowb = args.k;
    अगर (transb & 1) nrowb = args.n;

    info = -1;

    अगर (args.ldc < args.m) info = 13;
    अगर (args.ldb < nrowb)  info = 10;
    अगर (args.lda < nrowa)  info =  8;
    अगर (args.k < 0)        info =  5;
    अगर (args.n < 0)        info =  4;
    अगर (args.m < 0)        info =  3;
    अगर (transb < 0)        info =  2;
    अगर (transa < 0)        info =  1;
  }

  अगर (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (व्योम *)b;
    args.b = (व्योम *)a;
    args.c = (व्योम *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    अगर (TransB == CblasNoTrans)     transa = 0;
    अगर (TransB == CblasTrans)       transa = 1;
#अगर_अघोषित COMPLEX
    अगर (TransB == CblasConjNoTrans) transa = 0;
    अगर (TransB == CblasConjTrans)   transa = 1;
#अन्यथा
    अगर (TransB == CblasConjNoTrans) transa = 2;
    अगर (TransB == CblasConjTrans)   transa = 3;
#पूर्ण_अगर
    अगर (TransA == CblasNoTrans)     transb = 0;
    अगर (TransA == CblasTrans)       transb = 1;
#अगर_अघोषित COMPLEX
    अगर (TransA == CblasConjNoTrans) transb = 0;
    अगर (TransA == CblasConjTrans)   transb = 1;
#अन्यथा
    अगर (TransA == CblasConjNoTrans) transb = 2;
    अगर (TransA == CblasConjTrans)   transb = 3;
#पूर्ण_अगर

    nrowa = args.m;
    अगर (transa & 1) nrowa = args.k;
    nrowb = args.k;
    अगर (transb & 1) nrowb = args.n;

    info = -1;

    अगर (args.ldc < args.m) info = 13;
    अगर (args.ldb < nrowb)  info = 10;
    अगर (args.lda < nrowa)  info =  8;
    अगर (args.k < 0)        info =  5;
    अगर (args.n < 0)        info =  4;
    अगर (args.m < 0)        info =  3;
    अगर (transb < 0)        info =  2;
    अगर (transa < 0)        info =  1;

  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }


  अगर ((args.m == 0) || (args.n == 0)) वापस;
#अगर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#अगर defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
  अगर (support_avx512() )
#पूर्ण_अगर
  अगर (order == CblasRowMajor && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        वापस;
  }
#पूर्ण_अगर
#अगर defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
अगर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#अगर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#पूर्ण_अगर
)
#पूर्ण_अगर

  अगर (order == CblasRowMajor && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
        वापस;
  }
अन्यथा अगर (order == CblasRowMajor && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
        SGEMM_DIRECT_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        वापस;
  }

#पूर्ण_अगर //defined arm64
#पूर्ण_अगर //defined complex


#पूर्ण_अगर //defined CBLAS

#अगर !defined(BFLOAT16)  && !defined(HFLOAT16)
#अगर defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
अगर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#अगर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#पूर्ण_अगर
)
#पूर्ण_अगर //defined dynarch
{
अक्षर* TA,*TB;
      अगर (transa & 1) 
        TA = "T";
	अन्यथा 
	TA= "N";
      अगर (transb & 1) 
        TB = "T";
	अन्यथा 
	TB= "N";
#अगर_अघोषित COMPLEX
      अगर (transa == 3) 
	TA= "T";
      अगर (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#अगर_अघोषित DOUBLE
        SME_SGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#अन्यथा
        SME_DGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al, args.a, args.lda, args.b, args.ldb, be, args.c, args.ldc);
#पूर्ण_अगर
#अन्यथा
      अगर (transa == 2) 
	TA= "R";
      अगर (transb == 2) 
	TB= "R";
      अगर (transa == 3) 
	TA= "C";
      अगर (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)args.alpha;
FLOAT* be=(FLOAT*)args.beta;
#अगर_अघोषित DOUBLE
भग्न _Complex c_al={al[0],al[1]};
भग्न _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#अन्यथा
द्विगुन _Complex c_al={al[0],al[1]};
द्विगुन _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, args.m, args.n, args.k, al[0],al[1], args.a, args.lda, args.b, args.ldb, be[0],be[1], args.c, args.ldc);
#पूर्ण_अगर
#पूर्ण_अगर
	वापस;
}
#पूर्ण_अगर //defined arm64

#पूर्ण_अगर //defined b/hfloat16 

#अगर defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#अगर defined(DYNAMIC_ARCH)
  अगर (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    वापस;
  }
#पूर्ण_अगर
#अगर !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  अगर (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    वापस;
  }
#पूर्ण_अगर
#पूर्ण_अगर  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  अगर ((args.m == 0) || (args.n == 0)) वापस;

#अगर 0
  ख_लिखो(मानक_त्रुटि, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 args.m, args.n, args.k, args.lda, args.ldb, args.ldc);
#पूर्ण_अगर

#अगर (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#घोषणा BFLOAT16_GEMM_GEMV_FORWARD 1
#अन्यथा
#घोषणा BFLOAT16_GEMM_GEMV_FORWARD 0
#पूर्ण_अगर
#अगर (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#घोषणा HFLOAT16_GEMM_GEMV_FORWARD 1
#अन्यथा
#घोषणा HFLOAT16_GEMM_GEMV_FORWARD 0
#पूर्ण_अगर

#अगर defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#अगर defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#अन्यथा
  bool have_tuned_gemv = true;
#पूर्ण_अगर
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  अटल स्थिर अक्षर gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  अगर (args.k != 0) {
    अगर (args.n == 1) {
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
      अगर (transa & 1) {
        m = args.k;
        n = args.m;
      }
      अगर (transb & 1) {
        inc_x = args.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      अगर (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
        वापस;
      }
    }
    अगर (args.m == 1) {
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
      अगर (transa & 1) {
        inc_x = 1;
      }
      अगर (transb & 1) {
        m = args.n;
        n = args.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      अगर (is_efficient_gemv) {
        GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
        वापस;
      }
    }
  }
#पूर्ण_अगर

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#अगर USE_SMALL_MATRIX_OPT
अगर (*(FLOAT*)args.alpha != 0.0) {
#अगर !defined(COMPLEX)
  अगर(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){

	  अगर(*(FLOAT *)(args.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
	  }अन्यथा{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
	  }
	  वापस;
  }
#अन्यथा
  अगर(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, alpha[0], alpha[1], beta[0], beta[1])){
	  अगर(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, args.c, args.ldc);
	  }अन्यथा{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, alpha[0], alpha[1], args.b, args.ldb, beta[0], beta[1], args.c, args.ldc);
	  }
	  वापस;
  }
#पूर्ण_अगर
}
#पूर्ण_अगर

  बफर = (XFLOAT *)blas_memory_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#अगर defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)बफर + (WhereAmI() & 0xf) * GEMM_OFFSET_A);
#अन्यथा
  sa = (XFLOAT *)((BLASLONG)बफर +GEMM_OFFSET_A);
#पूर्ण_अगर
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#अगर_घोषित SMP
#अगर defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);
#पूर्ण_अगर

  MNK = (द्विगुन) args.m * (द्विगुन) args.n * (द्विगुन) args.k;
  args.nthreads = get_gemm_optimal_nthreads(MNK);

  args.common = शून्य;

 अगर (args.nthreads == 1) {
#पूर्ण_अगर

    (gemm[(transb << 2) | transa])(&args, शून्य, शून्य, sa, sb, 0);

#अगर_घोषित SMP

  } अन्यथा {

#अगर_अघोषित USE_SIMPLE_THREADED_LEVEL3

#अगर_अघोषित NO_AFFINITY
      nodes = get_अंक_nodes();

      अगर ((nodes > 1) && get_गाँठ_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, शून्य, शून्य, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } अन्यथा {
#पूर्ण_अगर

	(gemm[16 | (transb << 2) | transa])(&args, शून्य, शून्य, sa, sb, 0);

#अन्यथा

	GEMM_THREAD(mode, &args, शून्य, शून्य, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

#पूर्ण_अगर

#अगर_अघोषित USE_SIMPLE_THREADED_LEVEL3
#अगर_अघोषित NO_AFFINITY
      }
#पूर्ण_अगर
#पूर्ण_अगर

#पूर्ण_अगर

#अगर_घोषित SMP
  }
#पूर्ण_अगर

 blas_memory_free(बफर);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n, 2 * args.m * args.n * args.k);

  IDEBUG_END;

  वापस;
}
