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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ <stdlib.h>
#అన్తర్భావ <stdbool.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచించు SMP_పరిమితి_MIN 65536.0
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGEMM "
#నిర్వచించు GEMV BLASFUNC(qgemv)
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGEMM "
#నిర్వచించు GEMV BLASFUNC(dgemv)
#లేదావేళ defined(BFLOAT16)
#నిర్వచితమైతే BGEMM
#నిర్వచించు ERROR_NAME "BGEMM "
#నిర్వచించు GEMV BLASFUNC(bgemv)
#లేకపోతే
#నిర్వచించు ERROR_NAME "SBGEMM "
#నిర్వచించు GEMV BLASFUNC(sbgemv)
#అయితే_ముగింపు
#లేదావేళ defined(HFLOAT16)
#నిర్వచించు ERROR_NAME "SHGEMM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGEMM "
#నిర్వచించు GEMV BLASFUNC(sgemv)
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచించు SMP_పరిమితి_MIN 8192.0
#నిర్వచితం_కాకపోతే GEMM3M
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGEMM "
#నిర్వచించు GEMV BLASFUNC(xgemv)
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGEMM "
#నిర్వచించు GEMV BLASFUNC(zgemv)
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGEMM "
#నిర్వచించు GEMV BLASFUNC(cgemv)
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGEMM3M "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGEMM3M "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGEMM3M "
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే GEMM_MULTITHREAD_పరిమితి
#నిర్వచించు GEMM_MULTITHREAD_పరిమితి 4
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
#ఒకవేళ !defined(GEMM3M) || defined(GENERIC)
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_NN, GEMM_THREAD_TN, GEMM_THREAD_RN, GEMM_THREAD_CN,
  GEMM_THREAD_NT, GEMM_THREAD_TT, GEMM_THREAD_RT, GEMM_THREAD_CT,
  GEMM_THREAD_NR, GEMM_THREAD_TR, GEMM_THREAD_RR, GEMM_THREAD_CR,
  GEMM_THREAD_NC, GEMM_THREAD_TC, GEMM_THREAD_RC, GEMM_THREAD_CC,
#అయితే_ముగింపు
#లేకపోతే
  GEMM3M_NN, GEMM3M_TN, GEMM3M_RN, GEMM3M_CN,
  GEMM3M_NT, GEMM3M_TT, GEMM3M_RT, GEMM3M_CT,
  GEMM3M_NR, GEMM3M_TR, GEMM3M_RR, GEMM3M_CR,
  GEMM3M_NC, GEMM3M_TC, GEMM3M_RC, GEMM3M_CC,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM3M_THREAD_NN, GEMM3M_THREAD_TN, GEMM3M_THREAD_RN, GEMM3M_THREAD_CN,
  GEMM3M_THREAD_NT, GEMM3M_THREAD_TT, GEMM3M_THREAD_RT, GEMM3M_THREAD_CT,
  GEMM3M_THREAD_NR, GEMM3M_THREAD_TR, GEMM3M_THREAD_RR, GEMM3M_THREAD_CR,
  GEMM3M_THREAD_NC, GEMM3M_THREAD_TC, GEMM3M_THREAD_RC, GEMM3M_THREAD_CC,
#అయితే_ముగింపు
#అయితే_ముగింపు
};

#ఒకవేళ defined(SMALL_మాతృక_OPT) && !defined(GEMM3M) && !defined(XDOUBLE) && !defined(HFLOAT16) && !defined(BGEMM)
#నిర్వచించు USE_SMALL_మాతృక_OPT 1
#లేకపోతే
#నిర్వచించు USE_SMALL_మాతృక_OPT 0
#అయితే_ముగింపు

#ఒకవేళ USE_SMALL_మాతృక_OPT
#నిర్వచితం_కాకపోతే DYNAMIC_ARCH
#నిర్వచించు SMALL_KERNEL_ADDR(పట్టిక, idx) ((శూన్య *)(పట్టిక[idx]))
#లేకపోతే
#నిర్వచించు SMALL_KERNEL_ADDR(పట్టిక, idx) ((శూన్య *)(*(uintptr_t *)((అక్షరం *)gotoblas + (size_t)(పట్టిక[idx]))))
#అయితే_ముగింపు


#నిర్వచితం_కాకపోతే COMPLEX
స్థిర size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


స్థిర size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#నిర్వచించు GEMM_SMALL_KERNEL_B0(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#నిర్వచించు GEMM_SMALL_KERNEL(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#లేకపోతే

స్థిర size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

స్థిర size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#నిర్వచించు ZGEMM_SMALL_KERNEL(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#నిర్వచించు ZGEMM_SMALL_KERNEL_B0(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#అయితే_ముగింపు
#అయితే_ముగింపు

#ఒకవేళ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#నిర్వచించు XFEATURE_XTILEDATA 18
#నిర్వచించు ARCH_REQ_XCOMP_PERM 0x1023
స్థిర పూర్ణాంకం openblas_amxtile_permission = 0;
స్థిర పూర్ణాంకం init_amxtile_permission() {
  దీర్ఘ స్థితి_2 =
      syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  ఒకవేళ (స్థితి_2 != 0) {
    fprintf(stderr, "XTILEDATA permission not granted in your device(Linux, "
                    "Intel Sapphier Rapids), skip sbgemm calculation\n");
    తిరిగి -1;
  }
  openblas_amxtile_permission = 1;
  తిరిగి 0;
}
#అయితే_ముగింపు

#నిర్వచితమైతే DYNAMIC_ARCH
బాహ్య అక్షరం* gotoblas_corename(శూన్య);
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
#ఒకవేళ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
స్థిర అంతఃవరుస పూర్ణాంకం get_gemm_optimal_nthreads_neoversev1(ద్విగుణ MNK, పూర్ణాంకం ncpu) {
  తిరిగి
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
#అయితే_ముగింపు

#ఒకవేళ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
స్థిర అంతఃవరుస పూర్ణాంకం get_gemm_optimal_nthreads_neoversev2(ద్విగుణ MNK, పూర్ణాంకం ncpu) {
  తిరిగి
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
#అయితే_ముగింపు

స్థిర అంతఃవరుస పూర్ణాంకం get_gemm_optimal_nthreads(ద్విగుణ MNK) {
  పూర్ణాంకం ncpu = అంకె_cpu_avail(3);
#ఒకవేళ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  తిరిగి get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
#లేదావేళ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  తిరిగి get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
#లేదావేళ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16)
  ఒకవేళ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    తిరిగి get_gemm_optimal_nthreads_neoversev1(MNK, ncpu);
  }
  ఒకవేళ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    తిరిగి get_gemm_optimal_nthreads_neoversev2(MNK, ncpu);
  }
#అయితే_ముగింపు
  ఒకవేళ ( MNK <= (SMP_పరిమితి_MIN  * (ద్విగుణ) GEMM_MULTITHREAD_పరిమితి) ) {
    తిరిగి 1;
  }
  లేకపోతే {
    ఒకవేళ (MNK/ncpu < SMP_పరిమితి_MIN*(ద్విగుణ)GEMM_MULTITHREAD_పరిమితి) {
      తిరిగి MNK/(SMP_పరిమితి_MIN*(ద్విగుణ)GEMM_MULTITHREAD_పరిమితి);
    }
    లేకపోతే {
      తిరిగి ncpu;
    }
  }
}
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *TRANSA, అక్షరం *TRANSB,
	  blasint *M, blasint *N, blasint *K,
	  FLOAT *alpha,
	  IFLOAT *a, blasint *ldA,
	  IFLOAT *b, blasint *ldB,
	  FLOAT *beta,
	  FLOAT *c, blasint *ldC){

  blas_arg_t వాదనలు;

  పూర్ణాంకం transa, transb, nrowa, nrowb;
  blasint info;

  అక్షరం transA, transB;
  IFLOAT *బఫరు;
  IFLOAT *sa, *sb;

#నిర్వచితమైతే SMP
  ద్విగుణ MNK;
#ఒకవేళ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#ఒకవేళ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  పూర్ణాంకం nodes;
#అయితే_ముగింపు

  PRINT_DEBUG_NAME;

  వాదనలు.m = *M;
  వాదనలు.n = *N;
  వాదనలు.k = *K;

  వాదనలు.a = (శూన్య *)a;
  వాదనలు.b = (శూన్య *)b;
  వాదనలు.c = (శూన్య *)c;

  వాదనలు.lda = *ldA;
  వాదనలు.ldb = *ldB;
  వాదనలు.ldc = *ldC;

  వాదనలు.alpha = (శూన్య *)alpha;
  వాదనలు.beta  = (శూన్య *)beta;
  
  transA = *TRANSA;
  transB = *TRANSB;

  TOUPPER(transA);
  TOUPPER(transB);

  transa = -1;
  transb = -1;

  ఒకవేళ (transA == 'N') transa = 0;
  ఒకవేళ (transA == 'T') transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
  ఒకవేళ (transA == 'R') transa = 0;
  ఒకవేళ (transA == 'C') transa = 1;
#లేకపోతే
  ఒకవేళ (transA == 'R') transa = 2;
  ఒకవేళ (transA == 'C') transa = 3;
#అయితే_ముగింపు

  ఒకవేళ (transB == 'N') transb = 0;
  ఒకవేళ (transB == 'T') transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
  ఒకవేళ (transB == 'R') transb = 0;
  ఒకవేళ (transB == 'C') transb = 1;
#లేకపోతే
  ఒకవేళ (transB == 'R') transb = 2;
  ఒకవేళ (transB == 'C') transb = 3;
#అయితే_ముగింపు
  nrowa = వాదనలు.m;
  ఒకవేళ (transa & 1) nrowa = వాదనలు.k;
  nrowb = వాదనలు.k;
  ఒకవేళ (transb & 1) nrowb = వాదనలు.n;

  info = 0;

  ఒకవేళ (వాదనలు.ldc < వాదనలు.m) info = 13;
  ఒకవేళ (వాదనలు.ldb < nrowb)  info = 10;
  ఒకవేళ (వాదనలు.lda < nrowa)  info =  8;
  ఒకవేళ (వాదనలు.k < 0)        info =  5;
  ఒకవేళ (వాదనలు.n < 0)        info =  4;
  ఒకవేళ (వాదనలు.m < 0)        info =  3;
  ఒకవేళ (transb < 0)        info =  2;
  ఒకవేళ (transa < 0)        info =  1;

  ఒకవేళ (info){
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_TRANSPOSE TransA, గణన CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
#నిర్వచితం_కాకపోతే COMPLEX
	   FLOAT alpha,
	   IFLOAT *a, blasint lda,
	   IFLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#లేకపోతే
	   శూన్య *valpha,
	   శూన్య *va, blasint lda,
	   శూన్య *vb, blasint ldb,
	   శూన్య *vbeta,
	   శూన్య *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;
#అయితే_ముగింపు

  blas_arg_t వాదనలు;
  పూర్ణాంకం transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *బఫరు;
  XFLOAT *sa, *sb;

#నిర్వచితమైతే SMP
  ద్విగుణ MNK;
#ఒకవేళ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#ఒకవేళ defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  పూర్ణాంకం nodes;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

#నిర్వచితం_కాకపోతే COMPLEX
  వాదనలు.alpha = (శూన్య *)&alpha;
  వాదనలు.beta  = (శూన్య *)&beta;
#లేకపోతే
  వాదనలు.alpha = (శూన్య *)alpha;
  వాదనలు.beta  = (శూన్య *)beta;
#అయితే_ముగింపు

  transa = -1;
  transb = -1;
  info   =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    వాదనలు.m = m;
    వాదనలు.n = n;
    వాదనలు.k = k;

    వాదనలు.a = (శూన్య *)a;
    వాదనలు.b = (శూన్య *)b;
    వాదనలు.c = (శూన్య *)c;

    వాదనలు.lda = lda;
    వాదనలు.ldb = ldb;
    వాదనలు.ldc = ldc;

    ఒకవేళ (TransA == CblasNoTrans)     transa = 0;
    ఒకవేళ (TransA == CblasTrans)       transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (TransA == CblasConjNoTrans) transa = 0;
    ఒకవేళ (TransA == CblasConjTrans)   transa = 1;
#లేకపోతే
    ఒకవేళ (TransA == CblasConjNoTrans) transa = 2;
    ఒకవేళ (TransA == CblasConjTrans)   transa = 3;
#అయితే_ముగింపు
    ఒకవేళ (TransB == CblasNoTrans)     transb = 0;
    ఒకవేళ (TransB == CblasTrans)       transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (TransB == CblasConjNoTrans) transb = 0;
    ఒకవేళ (TransB == CblasConjTrans)   transb = 1;
#లేకపోతే
    ఒకవేళ (TransB == CblasConjNoTrans) transb = 2;
    ఒకవేళ (TransB == CblasConjTrans)   transb = 3;
#అయితే_ముగింపు

    nrowa = వాదనలు.m;
    ఒకవేళ (transa & 1) nrowa = వాదనలు.k;
    nrowb = వాదనలు.k;
    ఒకవేళ (transb & 1) nrowb = వాదనలు.n;

    info = -1;

    ఒకవేళ (వాదనలు.ldc < వాదనలు.m) info = 13;
    ఒకవేళ (వాదనలు.ldb < nrowb)  info = 10;
    ఒకవేళ (వాదనలు.lda < nrowa)  info =  8;
    ఒకవేళ (వాదనలు.k < 0)        info =  5;
    ఒకవేళ (వాదనలు.n < 0)        info =  4;
    ఒకవేళ (వాదనలు.m < 0)        info =  3;
    ఒకవేళ (transb < 0)        info =  2;
    ఒకవేళ (transa < 0)        info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    వాదనలు.m = n;
    వాదనలు.n = m;
    వాదనలు.k = k;

    వాదనలు.a = (శూన్య *)b;
    వాదనలు.b = (శూన్య *)a;
    వాదనలు.c = (శూన్య *)c;

    వాదనలు.lda = ldb;
    వాదనలు.ldb = lda;
    వాదనలు.ldc = ldc;

    ఒకవేళ (TransB == CblasNoTrans)     transa = 0;
    ఒకవేళ (TransB == CblasTrans)       transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (TransB == CblasConjNoTrans) transa = 0;
    ఒకవేళ (TransB == CblasConjTrans)   transa = 1;
#లేకపోతే
    ఒకవేళ (TransB == CblasConjNoTrans) transa = 2;
    ఒకవేళ (TransB == CblasConjTrans)   transa = 3;
#అయితే_ముగింపు
    ఒకవేళ (TransA == CblasNoTrans)     transb = 0;
    ఒకవేళ (TransA == CblasTrans)       transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (TransA == CblasConjNoTrans) transb = 0;
    ఒకవేళ (TransA == CblasConjTrans)   transb = 1;
#లేకపోతే
    ఒకవేళ (TransA == CblasConjNoTrans) transb = 2;
    ఒకవేళ (TransA == CblasConjTrans)   transb = 3;
#అయితే_ముగింపు

    nrowa = వాదనలు.m;
    ఒకవేళ (transa & 1) nrowa = వాదనలు.k;
    nrowb = వాదనలు.k;
    ఒకవేళ (transb & 1) nrowb = వాదనలు.n;

    info = -1;

    ఒకవేళ (వాదనలు.ldc < వాదనలు.m) info = 13;
    ఒకవేళ (వాదనలు.ldb < nrowb)  info = 10;
    ఒకవేళ (వాదనలు.lda < nrowa)  info =  8;
    ఒకవేళ (వాదనలు.k < 0)        info =  5;
    ఒకవేళ (వాదనలు.n < 0)        info =  4;
    ఒకవేళ (వాదనలు.m < 0)        info =  3;
    ఒకవేళ (transb < 0)        info =  2;
    ఒకవేళ (transa < 0)        info =  1;

  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }


  ఒకవేళ ((వాదనలు.m == 0) || (వాదనలు.n == 0)) తిరిగి;
#ఒకవేళ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ఒకవేళ defined(ARCH_x86) && (defined(USE_SGEMM_KERNEL_ప్రత్యక్ష)||defined(DYNAMIC_ARCH))
#ఒకవేళ defined(DYNAMIC_ARCH)
  ఒకవేళ (support_avx512() )
#అయితే_ముగింపు
  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ప్రత్యక్ష_PERFORMANT(m,n,k)) {
        SGEMM_ప్రత్యక్ష(m, n, k, a, lda, b, ldb, c, ldc);
        తిరిగి;
  }
#అయితే_ముగింపు
#ఒకవేళ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_ప్రత్యక్ష)||defined(DYNAMIC_ARCH))
#ఒకవేళ defined(DYNAMIC_ARCH)
ఒకవేళ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ఒకవేళ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#అయితే_ముగింపు
)
#అయితే_ముగింపు

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major && k==lda && n==ldb && n==ldc && beta == 0 && alpha == 1.0 && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ప్రత్యక్ష_PERFORMANT(m,n,k)) {
        SGEMM_ప్రత్యక్ష(m, n, k, a, lda, b, ldb, c, ldc);
        తిరిగి;
  }
లేకపోతే ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major && k==lda && n==ldb && n==ldc && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_ప్రత్యక్ష_PERFORMANT(m,n,k)) {
        SGEMM_ప్రత్యక్ష_ALPHA_BETA(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        తిరిగి;
  }

#అయితే_ముగింపు //defined arm64
#అయితే_ముగింపు //defined complex


#అయితే_ముగింపు //defined CBLAS

#ఒకవేళ !defined(BFLOAT16)  && !defined(HFLOAT16)
#ఒకవేళ defined(ARCH_ARM64) && (defined(USE_SGEMM_KERNEL_ప్రత్యక్ష)||defined(DYNAMIC_ARCH))
#ఒకవేళ defined(DYNAMIC_ARCH)
ఒకవేళ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ఒకవేళ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#అయితే_ముగింపు
)
#అయితే_ముగింపు //defined dynarch
{
అక్షరం* TA,*TB;
      ఒకవేళ (transa & 1) 
        TA = "T";
	లేకపోతే 
	TA= "N";
      ఒకవేళ (transb & 1) 
        TB = "T";
	లేకపోతే 
	TB= "N";
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (transa == 3) 
	TA= "T";
      ఒకవేళ (transb == 3) 
	TB= "T";
FLOAT* al=(FLOAT*)వాదనలు.alpha;
FLOAT* be=(FLOAT*)వాదనలు.beta;
#నిర్వచితం_కాకపోతే DOUBLE
        SME_SGEMM_KERNEL(TA,TB, వాదనలు.m, వాదనలు.n, వాదనలు.k, al, వాదనలు.a, వాదనలు.lda, వాదనలు.b, వాదనలు.ldb, be, వాదనలు.c, వాదనలు.ldc);
#లేకపోతే
        SME_DGEMM_KERNEL(TA,TB, వాదనలు.m, వాదనలు.n, వాదనలు.k, al, వాదనలు.a, వాదనలు.lda, వాదనలు.b, వాదనలు.ldb, be, వాదనలు.c, వాదనలు.ldc);
#అయితే_ముగింపు
#లేకపోతే
      ఒకవేళ (transa == 2) 
	TA= "R";
      ఒకవేళ (transb == 2) 
	TB= "R";
      ఒకవేళ (transa == 3) 
	TA= "C";
      ఒకవేళ (transb == 3) 
	TB= "C";
FLOAT* al=(FLOAT*)వాదనలు.alpha;
FLOAT* be=(FLOAT*)వాదనలు.beta;
#నిర్వచితం_కాకపోతే DOUBLE
దశాంశం _Complex c_al={al[0],al[1]};
దశాంశం _Complex c_be={be[0],be[1]};
        SME_CGEMM_KERNEL(TA,TB, వాదనలు.m, వాదనలు.n, వాదనలు.k, al[0],al[1], వాదనలు.a, వాదనలు.lda, వాదనలు.b, వాదనలు.ldb, be[0],be[1], వాదనలు.c, వాదనలు.ldc);
#లేకపోతే
ద్విగుణ _Complex c_al={al[0],al[1]};
ద్విగుణ _Complex c_be={be[0],be[1]};
        SME_ZGEMM_KERNEL(TA,TB, వాదనలు.m, వాదనలు.n, వాదనలు.k, al[0],al[1], వాదనలు.a, వాదనలు.lda, వాదనలు.b, వాదనలు.ldb, be[0],be[1], వాదనలు.c, వాదనలు.ldc);
#అయితే_ముగింపు
#అయితే_ముగింపు
	తిరిగి;
}
#అయితే_ముగింపు //defined arm64

#అయితే_ముగింపు //defined b/hfloat16 

#ఒకవేళ defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)
#ఒకవేళ defined(DYNAMIC_ARCH)
  ఒకవేళ (gotoblas->need_amxtile_permission &&
      openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    తిరిగి;
  }
#అయితే_ముగింపు
#ఒకవేళ !defined(DYNAMIC_ARCH) && defined(SAPPHIRERAPIDS)
  ఒకవేళ (openblas_amxtile_permission == 0 && init_amxtile_permission() == -1) {
    తిరిగి;
  }
#అయితే_ముగింపు
#అయితే_ముగింపు  // defined(__linux__) && defined(__x86_64__) && defined(BFLOAT16)

  ఒకవేళ ((వాదనలు.m == 0) || (వాదనలు.n == 0)) తిరిగి;

#ఒకవేళ 0
  fprintf(stderr, "m = %4d  n = %d  k = %d  lda = %4d  ldb = %4d  ldc = %4d\n",
	 వాదనలు.m, వాదనలు.n, వాదనలు.k, వాదనలు.lda, వాదనలు.ldb, వాదనలు.ldc);
#అయితే_ముగింపు

#ఒకవేళ (!defined(BFLOAT16) || (!defined(BGEMM) && defined(SBGEMM_GEMV_FORWARD)) || (defined(BGEMM) && defined(BGEMM_GEMV_FORWARD)))
#నిర్వచించు BFLOAT16_GEMM_GEMV_FORWARD 1
#లేకపోతే
#నిర్వచించు BFLOAT16_GEMM_GEMV_FORWARD 0
#అయితే_ముగింపు
#ఒకవేళ (!defined(HFLOAT16) || (!defined(HGEMM) && defined(SHGEMM_GEMV_FORWARD)) || (defined(HGEMM) && defined(HGEMM_GEMV_FORWARD)))
#నిర్వచించు HFLOAT16_GEMM_GEMV_FORWARD 1
#లేకపోతే
#నిర్వచించు HFLOAT16_GEMM_GEMV_FORWARD 0
#అయితే_ముగింపు

#ఒకవేళ defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
#ఒకవేళ defined(ARCH_ARM64)
  // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
  // perform poorly in certain circumstances. We use the following boolean
  // variable along with the gemv argument values to avoid these inefficient
  // gemv cases, see github issue#4951.
  bool have_tuned_gemv = false;
#లేకపోతే
  bool have_tuned_gemv = true;
#అయితే_ముగింపు
  // GEMV transpose modes indexed by the GEMM transpose codes (N=0, T=1, R=2, C=3).
  // For real types only N and T occur. For complex types the OpenBLAS GEMV
  // extensions O, U, S and D (index 4..7) behave like N, T, R and C but
  // additionally conjugate the vector x, which lets a conjugated vector-shaped
  // operand of GEMM be handled without an explicit conjugation pass.
  // Even indices run on the GEMV_N kernels, odd indices on the GEMV_T kernels.
  స్థిర స్థిరరాశి అక్షరం gemv_trans[] = "NTRCOUSD";
  // Check if we can convert GEMM -> GEMV
  ఒకవేళ (వాదనలు.k != 0) {
    ఒకవేళ (వాదనలు.n == 1) {
      blasint inc_x = 1;
      blasint inc_y = 1;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = వాదనలు.m;
      blasint n = వాదనలు.k;
      blasint lda = వాదనలు.lda;
      // C = alpha * op(A) * op(B) + beta * C with op(B) a k x 1 vector: GEMV on A,
      // where a conjugated op(B) selects the x-conjugating modes.
      పూర్ణాంకం trans_idx = transa | ((transb & 2) << 1);
      అక్షరం NT = gemv_trans[trans_idx];
      ఒకవేళ (transa & 1) {
        m = వాదనలు.k;
        n = వాదనలు.m;
      }
      ఒకవేళ (transb & 1) {
        inc_x = వాదనలు.ldb;
      }
      bool is_efficient_gemv = have_tuned_gemv || (!(trans_idx & 1) || inc_x == 1);
      ఒకవేళ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, వాదనలు.alpha, వాదనలు.a, &lda, వాదనలు.b, &inc_x, వాదనలు.beta, వాదనలు.c, &inc_y);
        తిరిగి;
      }
    }
    ఒకవేళ (వాదనలు.m == 1) {
      blasint inc_x = వాదనలు.lda;
      blasint inc_y = వాదనలు.ldc;
      // These were passed in as blasint, but the struct translates them to blaslong
      blasint m = వాదనలు.k;
      blasint n = వాదనలు.n;
      blasint ldb = వాదనలు.ldb;
      // C^T = op(B)^T * op(A)^T with op(A) a 1 x k vector: GEMV on B, where the
      // transposition swaps N<->T and R<->C, and a conjugated op(A) selects the
      // x-conjugating modes.
      పూర్ణాంకం trans_idx = (transb ^ 1) | ((transa & 2) << 1);
      అక్షరం NT = gemv_trans[trans_idx];
      ఒకవేళ (transa & 1) {
        inc_x = 1;
      }
      ఒకవేళ (transb & 1) {
        m = వాదనలు.n;
        n = వాదనలు.k;
      }
      bool is_efficient_gemv = have_tuned_gemv || ((!(trans_idx & 1) && inc_y == 1) || ((trans_idx & 1) && inc_x == 1));
      ఒకవేళ (is_efficient_gemv) {
        GEMV(&NT, &m, &n, వాదనలు.alpha, వాదనలు.b, &ldb, వాదనలు.a, &inc_x, వాదనలు.beta, వాదనలు.c, &inc_y);
        తిరిగి;
      }
    }
  }
#అయితే_ముగింపు

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

#ఒకవేళ USE_SMALL_మాతృక_OPT
ఒకవేళ (*(FLOAT*)వాదనలు.alpha != 0.0) {
#ఒకవేళ !defined(COMPLEX)
  ఒకవేళ(GEMM_SMALL_మాతృక_PERMIT(transa, transb, వాదనలు.m, వాదనలు.n, వాదనలు.k, *(FLOAT *)(వాదనలు.alpha), *(FLOAT *)(వాదనలు.beta))){

	  ఒకవేళ(*(FLOAT *)(వాదనలు.beta) == 0.0){
		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(వాదనలు.m, వాదనలు.n, వాదనలు.k, వాదనలు.a, వాదనలు.lda, *(FLOAT *)(వాదనలు.alpha), వాదనలు.b, వాదనలు.ldb, వాదనలు.c, వాదనలు.ldc);
	  }లేకపోతే{
		(GEMM_SMALL_KERNEL((transb << 2) | transa))(వాదనలు.m, వాదనలు.n, వాదనలు.k, వాదనలు.a, వాదనలు.lda, *(FLOAT *)(వాదనలు.alpha), వాదనలు.b, వాదనలు.ldb, *(FLOAT *)(వాదనలు.beta), వాదనలు.c, వాదనలు.ldc);
	  }
	  తిరిగి;
  }
#లేకపోతే
  ఒకవేళ(GEMM_SMALL_మాతృక_PERMIT(transa, transb, వాదనలు.m, వాదనలు.n, వాదనలు.k, alpha[0], alpha[1], beta[0], beta[1])){
	  ఒకవేళ(beta[0] == 0.0 && beta[1] == 0.0){
		(ZGEMM_SMALL_KERNEL_B0((transb << 2) | transa))(వాదనలు.m, వాదనలు.n, వాదనలు.k, వాదనలు.a, వాదనలు.lda, alpha[0], alpha[1], వాదనలు.b, వాదనలు.ldb, వాదనలు.c, వాదనలు.ldc);
	  }లేకపోతే{
		(ZGEMM_SMALL_KERNEL((transb << 2) | transa))(వాదనలు.m, వాదనలు.n, వాదనలు.k, వాదనలు.a, వాదనలు.lda, alpha[0], alpha[1], వాదనలు.b, వాదనలు.ldb, beta[0], beta[1], వాదనలు.c, వాదనలు.ldc);
	  }
	  తిరిగి;
  }
#అయితే_ముగింపు
}
#అయితే_ముగింపు

  బఫరు = (XFLOAT *)blas_జ్ఞాపకం_alloc(0);
//For LOONGARCH64, applying an offset to the buffer is essential
//for minimizing cache conflicts and optimizing performance.
#ఒకవేళ defined(ARCH_LOONGARCH64) && !defined(NO_AFFINITY)
  sa = (XFLOAT *)((BLASLONG)బఫరు + (WhereAmI() & 0xf) * GEMM_స్థానభ్రంశం_A);
#లేకపోతే
  sa = (XFLOAT *)((BLASLONG)బఫరు +GEMM_స్థానభ్రంశం_A);
#అయితే_ముగింపు
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);

#నిర్వచితమైతే SMP
#ఒకవేళ defined(USE_SIMPLE_THREADED_LEVEL3) || !defined(NO_AFFINITY)
  బహుళకం |= (transa << BLAS_TRANSA_జరుపు);
  బహుళకం |= (transb << BLAS_TRANSB_జరుపు);
#అయితే_ముగింపు

  MNK = (ద్విగుణ) వాదనలు.m * (ద్విగుణ) వాదనలు.n * (ద్విగుణ) వాదనలు.k;
  వాదనలు.nthreads = get_gemm_optimal_nthreads(MNK);

  వాదనలు.common = ఖాళీ;

 ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

    (gemm[(transb << 2) | transa])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#నిర్వచితమైతే SMP

  } లేకపోతే {

#నిర్వచితం_కాకపోతే USE_SIMPLE_THREADED_LEVEL3

#నిర్వచితం_కాకపోతే NO_AFFINITY
      nodes = get_అంకె_nodes();

      ఒకవేళ ((nodes > 1) && get_కణుపు_equal()) {

	వాదనలు.nthreads /= nodes;

	gemm_thread_mn(బహుళకం, &వాదనలు, ఖాళీ, ఖాళీ, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } లేకపోతే {
#అయితే_ముగింపు

	(gemm[16 | (transb << 2) | transa])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#లేకపోతే

	GEMM_THREAD(బహుళకం, &వాదనలు, ఖాళీ, ఖాళీ, gemm[(transb << 2) | transa], sa, sb, వాదనలు.nthreads);

#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే USE_SIMPLE_THREADED_LEVEL3
#నిర్వచితం_కాకపోతే NO_AFFINITY
      }
#అయితే_ముగింపు
#అయితే_ముగింపు

#అయితే_ముగింపు

#నిర్వచితమైతే SMP
  }
#అయితే_ముగింపు

 blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE, వాదనలు.m * వాదనలు.k + వాదనలు.k * వాదనలు.n + వాదనలు.m * వాదనలు.n, 2 * వాదనలు.m * వాదనలు.n * వాదనలు.k);

  IDEBUG_చివర;

  తిరిగి;
}
