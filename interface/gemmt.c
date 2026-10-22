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

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#وگرنهاگر defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#وگرنه
#define ERROR_NAME "SGEMMTR"
#endif
#وگرنه
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#وگرنهاگر defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#وگرنه
#define ERROR_NAME "SGEMMT "
#endif
#endif
#وگرنه
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#وگرنه
#define ERROR_NAME "CGEMMTR"
#endif
#وگرنه
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#وگرنه
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *UPLO, نویسه *TRANSA, نویسه *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	صحیح transa, transb, uplo;
	blasint info;

	نویسه transA, transB, Uplo;
	blasint nrowa, nrowb;
#اگر !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#اگر defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#وگرنه
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

	اگر (transA == 'N')
		transa = 0;
	اگر (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	اگر (transA == 'R')
		transa = 0;
	اگر (transA == 'C')
		transa = 1;
#وگرنه
	اگر (transA == 'R')
		transa = 2;
	اگر (transA == 'C')
		transa = 3;
#endif

	اگر (transB == 'N')
		transb = 0;
	اگر (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	اگر (transB == 'R')
		transb = 0;
	اگر (transB == 'C')
		transb = 1;
#وگرنه
	اگر (transB == 'R')
		transb = 2;
	اگر (transB == 'C')
		transb = 3;
#endif

	اگر (Uplo == 'U')
		uplo = 0;
	اگر (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	اگر (transa & 1) nrowa = k;
	nrowb = k;
	اگر (transb & 1) nrowb = m;

	info = 0;

	اگر (ldc < MAX(1, m))
		info = 13;
	اگر (ldb < MAX(1, nrowb))
		info = 10;
	اگر (lda < MAX(1, nrowa))
		info = 8;
	اگر (k < 0)
		info = 5;
	اگر (m < 0)
		info = 4;
	اگر (transb < 0)
		info = 3;
	اگر (transa < 0)
		info = 2;
	اگر (uplo < 0)
		info = 1;

	اگر (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, اندازه(ERROR_NAME));
		بازگشت;
	}
#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order, شمارش CBLAS_UPLO Uplo,
	   شمارش CBLAS_TRANSPOSE TransA, شمارش CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#وگرنه
	   تهی *valpha,
	   تهی *va, blasint LDA,
	   تهی *vb, blasint LDB, تهی *vbeta, تهی *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	صحیح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	اگر (order == CblasColMajor) {
		اگر (Uplo == CblasUpper) uplo = 0;
		اگر (Uplo == CblasLower) uplo = 1;

		اگر (TransA == CblasNoTrans)
			transa = 0;
		اگر (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		اگر (TransA == CblasConjNoTrans)
			transa = 0;
		اگر (TransA == CblasConjTrans)
			transa = 1;
#وگرنه
		اگر (TransA == CblasConjNoTrans)
			transa = 2;
		اگر (TransA == CblasConjTrans)
			transa = 3;
#endif
		اگر (TransB == CblasNoTrans)
			transb = 0;
		اگر (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		اگر (TransB == CblasConjNoTrans)
			transb = 0;
		اگر (TransB == CblasConjTrans)
			transb = 1;
#وگرنه
		اگر (TransB == CblasConjNoTrans)
			transb = 2;
		اگر (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (تهی *)A;
		b = (تهی *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		اگر (transa & 1) nrowa = k;
		nrowb = k;
		اگر (transb & 1) nrowb = m;

		اگر (ldc < MAX(1, m))
			info = 13;
		اگر (ldb < MAX(1, nrowb))
			info = 10;
		اگر (lda < MAX(1, nrowa))
			info = 8;
		اگر (k < 0)
			info = 5;
		اگر (m < 0)
			info = 4;
		اگر (transb < 0)
			info = 3;
		اگر (transa < 0)
			info = 2;
		اگر (uplo < 0)
			info = 1;
	}

	اگر (order == CblasRowMajor) {

		a = (تهی *)B;
		b = (تهی *)A;

		lda = LDB;
		ldb = LDA;

		اگر (Uplo == CblasUpper) uplo = 1;
		اگر (Uplo == CblasLower) uplo = 0;

		اگر (TransB == CblasNoTrans)
			transa = 0;
		اگر (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		اگر (TransB == CblasConjNoTrans)
			transa = 0;
		اگر (TransB == CblasConjTrans)
			transa = 1;
#وگرنه
		اگر (TransB == CblasConjNoTrans)
			transa = 2;
		اگر (TransB == CblasConjTrans)
			transa = 3;
#endif
		اگر (TransA == CblasNoTrans)
			transb = 0;
		اگر (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		اگر (TransA == CblasConjNoTrans)
			transb = 0;
		اگر (TransA == CblasConjTrans)
			transb = 1;
#وگرنه
		اگر (TransA == CblasConjNoTrans)
			transb = 2;
		اگر (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		اگر (transa & 1) ncola = k;
		ncolb = k;
		اگر (transb & 1) ncolb = m;

		اگر (ldc < MAX(1,m))
			info = 13;
		اگر (ldb < MAX(1, ncolb))
			info = 8;
		اگر (lda < MAX(1, ncola))
			info = 10;
		اگر (k < 0)
			info = 5;
		اگر (m < 0)
			info = 4;
		اگر (transb < 0)
			info = 2;
		اگر (transa < 0)
			info = 3;
		اگر (uplo < 0)
			info = 1;
	}

	اگر (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, اندازه(ERROR_NAME));
		بازگشت;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#وگرنه
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	دوگانه MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	صحیح mode = BLAS_XDOUBLE | BLAS_REAL;
#وگرنهاگر defined(DOUBLE)
	صحیح mode = BLAS_DOUBLE  | BLAS_REAL;
#وگرنه
	صحیح mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#وگرنه
#ifdef XDOUBLE
	صحیح mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنهاگر defined(DOUBLE)
	صحیح mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#وگرنه
	صحیح mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ایستا صحیح (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#وگرنه
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

	اگر (m == 0)
		بازگشت;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (تهی *)a;
	args.b = (تهی *)b;
	args.c = (تهی *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#اگر defined(COMPLEX)
	args.alpha = (تهی *)alpha;
	args.beta  = (تهی *)beta;
#وگرنه
	args.alpha = (تهی *)&alpha;
	args.beta  = (تهی *)&beta;
#endif

#ifdef SMP
	args.common   = پوچ;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	اگر (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, اندازه(ERROR_NAME));
		بازگشت;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (دوگانه)(m + 1) * (دوگانه)m * (دوگانه)k;
	اگر (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} وگرنه {
		args.nthreads = رقم_cpu_avail(3);
	}

	اگر (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, پوچ, پوچ, sa, sb, 0);

#ifdef SMP
	} وگرنه {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, پوچ, پوچ, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	بازگشت;
}
