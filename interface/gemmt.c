/*********************************************************************/
/* Copyright 2022, The OpenBLAS Project.                             */
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
/*********************************************************************/

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <stdlib.h>
#অন্তর্ভাব "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#নতুবা_যদি defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#অন্যথা
#define ERROR_NAME "SGEMMTR"
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#নতুবা_যদি defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#অন্যথা
#define ERROR_NAME "SGEMMT "
#endif
#endif
#অন্যথা
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#অন্যথা
#define ERROR_NAME "CGEMMTR"
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#অন্যথা
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(অক্ষর *UPLO, অক্ষর *TRANSA, অক্ষর *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	পূর্ণসংখ্যা transa, transb, uplo;
	blasint info;

	অক্ষর transA, transB, Uplo;
	blasint nrowa, nrowb;
#যদি !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#যদি defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#অন্যথা
	alpha = *Alpha;
	beta = *Beta;
#endif

	lda = *ldA;
	ldb = *ldB;
	ldc = *ldC;

	transA = *TRANSA;
	transB = *TRANSB;
	Uplo = *UPLO;
	TOUPPER(transA);
	TOUPPER(transB);
	TOUPPER(Uplo);

	transa = -1;
	transb = -1;
	uplo = -1;

	যদি (transA == 'N')
		transa = 0;
	যদি (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	যদি (transA == 'R')
		transa = 0;
	যদি (transA == 'C')
		transa = 1;
#অন্যথা
	যদি (transA == 'R')
		transa = 2;
	যদি (transA == 'C')
		transa = 3;
#endif

	যদি (transB == 'N')
		transb = 0;
	যদি (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	যদি (transB == 'R')
		transb = 0;
	যদি (transB == 'C')
		transb = 1;
#অন্যথা
	যদি (transB == 'R')
		transb = 2;
	যদি (transB == 'C')
		transb = 3;
#endif

	যদি (Uplo == 'U')
		uplo = 0;
	যদি (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	যদি (transa & 1) nrowa = k;
	nrowb = k;
	যদি (transb & 1) nrowb = m;

	info = 0;

	যদি (ldc < MAX(1, m))
		info = 13;
	যদি (ldb < MAX(1, nrowb))
		info = 10;
	যদি (lda < MAX(1, nrowa))
		info = 8;
	যদি (k < 0)
		info = 5;
	যদি (m < 0)
		info = 4;
	যদি (transb < 0)
		info = 3;
	যদি (transa < 0)
		info = 2;
	যদি (uplo < 0)
		info = 1;

	যদি (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, আকার(ERROR_NAME));
		ফেরত;
	}
#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo,
	   গণনা CBLAS_TRANSPOSE TransA, গণনা CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#অন্যথা
	   শূন্য *valpha,
	   শূন্য *va, blasint LDA,
	   শূন্য *vb, blasint LDB, শূন্য *vbeta, শূন্য *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	পূর্ণসংখ্যা transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	যদি (order == CblasColMajor) {
		যদি (Uplo == CblasUpper) uplo = 0;
		যদি (Uplo == CblasLower) uplo = 1;

		যদি (TransA == CblasNoTrans)
			transa = 0;
		যদি (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		যদি (TransA == CblasConjNoTrans)
			transa = 0;
		যদি (TransA == CblasConjTrans)
			transa = 1;
#অন্যথা
		যদি (TransA == CblasConjNoTrans)
			transa = 2;
		যদি (TransA == CblasConjTrans)
			transa = 3;
#endif
		যদি (TransB == CblasNoTrans)
			transb = 0;
		যদি (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		যদি (TransB == CblasConjNoTrans)
			transb = 0;
		যদি (TransB == CblasConjTrans)
			transb = 1;
#অন্যথা
		যদি (TransB == CblasConjNoTrans)
			transb = 2;
		যদি (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (শূন্য *)A;
		b = (শূন্য *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		যদি (transa & 1) nrowa = k;
		nrowb = k;
		যদি (transb & 1) nrowb = m;

		যদি (ldc < MAX(1, m))
			info = 13;
		যদি (ldb < MAX(1, nrowb))
			info = 10;
		যদি (lda < MAX(1, nrowa))
			info = 8;
		যদি (k < 0)
			info = 5;
		যদি (m < 0)
			info = 4;
		যদি (transb < 0)
			info = 3;
		যদি (transa < 0)
			info = 2;
		যদি (uplo < 0)
			info = 1;
	}

	যদি (order == CblasRowMajor) {

		a = (শূন্য *)B;
		b = (শূন্য *)A;

		lda = LDB;
		ldb = LDA;

		যদি (Uplo == CblasUpper) uplo = 1;
		যদি (Uplo == CblasLower) uplo = 0;

		যদি (TransB == CblasNoTrans)
			transa = 0;
		যদি (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		যদি (TransB == CblasConjNoTrans)
			transa = 0;
		যদি (TransB == CblasConjTrans)
			transa = 1;
#অন্যথা
		যদি (TransB == CblasConjNoTrans)
			transa = 2;
		যদি (TransB == CblasConjTrans)
			transa = 3;
#endif
		যদি (TransA == CblasNoTrans)
			transb = 0;
		যদি (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		যদি (TransA == CblasConjNoTrans)
			transb = 0;
		যদি (TransA == CblasConjTrans)
			transb = 1;
#অন্যথা
		যদি (TransA == CblasConjNoTrans)
			transb = 2;
		যদি (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		যদি (transa & 1) ncola = k;
		ncolb = k;
		যদি (transb & 1) ncolb = m;

		যদি (ldc < MAX(1,m))
			info = 13;
		যদি (ldb < MAX(1, ncolb))
			info = 8;
		যদি (lda < MAX(1, ncola))
			info = 10;
		যদি (k < 0)
			info = 5;
		যদি (m < 0)
			info = 4;
		যদি (transb < 0)
			info = 2;
		যদি (transa < 0)
			info = 3;
		যদি (uplo < 0)
			info = 1;
	}

	যদি (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, আকার(ERROR_NAME));
		ফেরত;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_সূচক ((uplo << 2) | (transb << 1) | transa)
#অন্যথা
#define GEMMT_TABLE_সূচক ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *বাফার, *sa, *sb;

#ifdef SMP
	দ্বিগুণ MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	পূর্ণসংখ্যা mode = BLAS_XDOUBLE | BLAS_REAL;
#নতুবা_যদি defined(DOUBLE)
	পূর্ণসংখ্যা mode = BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
	পূর্ণসংখ্যা mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
	পূর্ণসংখ্যা mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#নতুবা_যদি defined(DOUBLE)
	পূর্ণসংখ্যা mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
	পূর্ণসংখ্যা mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	স্থির পূর্ণসংখ্যা (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#অন্যথা
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#endif
	};

	যদি (m == 0)
		ফেরত;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (শূন্য *)a;
	args.b = (শূন্য *)b;
	args.c = (শূন্য *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#যদি defined(COMPLEX)
	args.alpha = (শূন্য *)alpha;
	args.beta  = (শূন্য *)beta;
#অন্যথা
	args.alpha = (শূন্য *)&alpha;
	args.beta  = (শূন্য *)&beta;
#endif

#ifdef SMP
	args.common   = ফাঁকা;
	args.nthreads = 1;
#endif

	বাফার = (FLOAT *)blas_memory_alloc(0);
	যদি (!বাফার) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, আকার(ERROR_NAME));
		ফেরত;
	}

	sa = (FLOAT *)((BLASLONG)বাফার + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * আকার_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (দ্বিগুণ)(m + 1) * (দ্বিগুণ)m * (দ্বিগুণ)k;
	যদি (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} অন্যথা {
		args.nthreads = অঙ্ক_cpu_avail(3);
	}

	যদি (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_সূচক])(&args, ফাঁকা, ফাঁকা, sa, sb, 0);

#ifdef SMP
	} অন্যথা {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, ফাঁকা, ফাঁকা, gemmt[GEMMT_TABLE_সূচক], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(বাফার);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	ফেরত;
}
