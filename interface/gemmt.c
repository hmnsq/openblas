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

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#وإلاإذا defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#ܐܚܪܢܐ
#define ERROR_NAME "SGEMMTR"
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#وإلاإذا defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#ܐܚܪܢܐ
#define ERROR_NAME "SGEMMT "
#endif
#endif
#ܐܚܪܢܐ
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#ܐܚܪܢܐ
#define ERROR_NAME "CGEMMTR"
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#ܐܚܪܢܐ
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANSA, حرف *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	صحيح transa, transb, uplo;
	blasint info;

	حرف transA, transB, Uplo;
	blasint nrowa, nrowb;
#ܐܢ !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ܐܢ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ܐܚܪܢܐ
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

	ܐܢ (transA == 'N')
		transa = 0;
	ܐܢ (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	ܐܢ (transA == 'R')
		transa = 0;
	ܐܢ (transA == 'C')
		transa = 1;
#ܐܚܪܢܐ
	ܐܢ (transA == 'R')
		transa = 2;
	ܐܢ (transA == 'C')
		transa = 3;
#endif

	ܐܢ (transB == 'N')
		transb = 0;
	ܐܢ (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	ܐܢ (transB == 'R')
		transb = 0;
	ܐܢ (transB == 'C')
		transb = 1;
#ܐܚܪܢܐ
	ܐܢ (transB == 'R')
		transb = 2;
	ܐܢ (transB == 'C')
		transb = 3;
#endif

	ܐܢ (Uplo == 'U')
		uplo = 0;
	ܐܢ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ܐܢ (transa & 1) nrowa = k;
	nrowb = k;
	ܐܢ (transb & 1) nrowb = m;

	info = 0;

	ܐܢ (ldc < MAX(1, m))
		info = 13;
	ܐܢ (ldb < MAX(1, nrowb))
		info = 10;
	ܐܢ (lda < MAX(1, nrowa))
		info = 8;
	ܐܢ (k < 0)
		info = 5;
	ܐܢ (m < 0)
		info = 4;
	ܐܢ (transb < 0)
		info = 3;
	ܐܢ (transa < 0)
		info = 2;
	ܐܢ (uplo < 0)
		info = 1;

	ܐܢ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		ܦܢܐ;
	}
#ܐܚܪܢܐ

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ܐܚܪܢܐ
	   فراغ *valpha,
	   فراغ *va, blasint LDA,
	   فراغ *vb, blasint LDB, فراغ *vbeta, فراغ *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	صحيح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ܐܢ (order == CblasColMajor) {
		ܐܢ (Uplo == CblasUpper) uplo = 0;
		ܐܢ (Uplo == CblasLower) uplo = 1;

		ܐܢ (TransA == CblasNoTrans)
			transa = 0;
		ܐܢ (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ܐܢ (TransA == CblasConjNoTrans)
			transa = 0;
		ܐܢ (TransA == CblasConjTrans)
			transa = 1;
#ܐܚܪܢܐ
		ܐܢ (TransA == CblasConjNoTrans)
			transa = 2;
		ܐܢ (TransA == CblasConjTrans)
			transa = 3;
#endif
		ܐܢ (TransB == CblasNoTrans)
			transb = 0;
		ܐܢ (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ܐܢ (TransB == CblasConjNoTrans)
			transb = 0;
		ܐܢ (TransB == CblasConjTrans)
			transb = 1;
#ܐܚܪܢܐ
		ܐܢ (TransB == CblasConjNoTrans)
			transb = 2;
		ܐܢ (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (فراغ *)A;
		b = (فراغ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ܐܢ (transa & 1) nrowa = k;
		nrowb = k;
		ܐܢ (transb & 1) nrowb = m;

		ܐܢ (ldc < MAX(1, m))
			info = 13;
		ܐܢ (ldb < MAX(1, nrowb))
			info = 10;
		ܐܢ (lda < MAX(1, nrowa))
			info = 8;
		ܐܢ (k < 0)
			info = 5;
		ܐܢ (m < 0)
			info = 4;
		ܐܢ (transb < 0)
			info = 3;
		ܐܢ (transa < 0)
			info = 2;
		ܐܢ (uplo < 0)
			info = 1;
	}

	ܐܢ (order == CblasRowMajor) {

		a = (فراغ *)B;
		b = (فراغ *)A;

		lda = LDB;
		ldb = LDA;

		ܐܢ (Uplo == CblasUpper) uplo = 1;
		ܐܢ (Uplo == CblasLower) uplo = 0;

		ܐܢ (TransB == CblasNoTrans)
			transa = 0;
		ܐܢ (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ܐܢ (TransB == CblasConjNoTrans)
			transa = 0;
		ܐܢ (TransB == CblasConjTrans)
			transa = 1;
#ܐܚܪܢܐ
		ܐܢ (TransB == CblasConjNoTrans)
			transa = 2;
		ܐܢ (TransB == CblasConjTrans)
			transa = 3;
#endif
		ܐܢ (TransA == CblasNoTrans)
			transb = 0;
		ܐܢ (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ܐܢ (TransA == CblasConjNoTrans)
			transb = 0;
		ܐܢ (TransA == CblasConjTrans)
			transb = 1;
#ܐܚܪܢܐ
		ܐܢ (TransA == CblasConjNoTrans)
			transb = 2;
		ܐܢ (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ܐܢ (transa & 1) ncola = k;
		ncolb = k;
		ܐܢ (transb & 1) ncolb = m;

		ܐܢ (ldc < MAX(1,m))
			info = 13;
		ܐܢ (ldb < MAX(1, ncolb))
			info = 8;
		ܐܢ (lda < MAX(1, ncola))
			info = 10;
		ܐܢ (k < 0)
			info = 5;
		ܐܢ (m < 0)
			info = 4;
		ܐܢ (transb < 0)
			info = 2;
		ܐܢ (transa < 0)
			info = 3;
		ܐܢ (uplo < 0)
			info = 1;
	}

	ܐܢ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		ܦܢܐ;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#ܐܚܪܢܐ
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	مزدوج MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	صحيح mode = BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
	صحيح mode = BLAS_DOUBLE  | BLAS_REAL;
#ܐܚܪܢܐ
	صحيح mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
	صحيح mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
	صحيح mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#ܐܚܪܢܐ
	صحيح mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ساكن صحيح (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ܐܚܪܢܐ
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

	ܐܢ (m == 0)
		ܦܢܐ;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (فراغ *)a;
	args.b = (فراغ *)b;
	args.c = (فراغ *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#ܐܢ defined(COMPLEX)
	args.alpha = (فراغ *)alpha;
	args.beta  = (فراغ *)beta;
#ܐܚܪܢܐ
	args.alpha = (فراغ *)&alpha;
	args.beta  = (فراغ *)&beta;
#endif

#ifdef SMP
	args.common   = عدم;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	ܐܢ (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		ܦܢܐ;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (مزدوج)(m + 1) * (مزدوج)m * (مزدوج)k;
	ܐܢ (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} ܐܚܪܢܐ {
		args.nthreads = num_cpu_avail(3);
	}

	ܐܢ (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP
	} ܐܚܪܢܐ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, عدم, عدم, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	ܦܢܐ;
}
