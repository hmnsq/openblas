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

#include <stdio.h>
#include <stdlib.h>
#include "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#elif defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#បើមិនដូច្នេះ
#define ERROR_NAME "SGEMMTR"
#endif
#បើមិនដូច្នេះ
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#elif defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#បើមិនដូច្នេះ
#define ERROR_NAME "SGEMMT "
#endif
#endif
#បើមិនដូច្នេះ
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#បើមិនដូច្នេះ
#define ERROR_NAME "CGEMMTR"
#endif
#បើមិនដូច្នេះ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#បើមិនដូច្នេះ
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANSA, var_na *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ចំនួនគត់ transa, transb, uplo;
	blasint info;

	var_na transA, transB, Uplo;
	blasint nrowa, nrowb;
#ប្រសិនបើ !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ប្រសិនបើ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#បើមិនដូច្នេះ
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

	ប្រសិនបើ (transA == 'N')
		transa = 0;
	ប្រសិនបើ (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	ប្រសិនបើ (transA == 'R')
		transa = 0;
	ប្រសិនបើ (transA == 'C')
		transa = 1;
#បើមិនដូច្នេះ
	ប្រសិនបើ (transA == 'R')
		transa = 2;
	ប្រសិនបើ (transA == 'C')
		transa = 3;
#endif

	ប្រសិនបើ (transB == 'N')
		transb = 0;
	ប្រសិនបើ (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	ប្រសិនបើ (transB == 'R')
		transb = 0;
	ប្រសិនបើ (transB == 'C')
		transb = 1;
#បើមិនដូច្នេះ
	ប្រសិនបើ (transB == 'R')
		transb = 2;
	ប្រសិនបើ (transB == 'C')
		transb = 3;
#endif

	ប្រសិនបើ (Uplo == 'U')
		uplo = 0;
	ប្រសិនបើ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ប្រសិនបើ (transa & 1) nrowa = k;
	nrowb = k;
	ប្រសិនបើ (transb & 1) nrowb = m;

	info = 0;

	ប្រសិនបើ (ldc < MAX(1, m))
		info = 13;
	ប្រសិនបើ (ldb < MAX(1, nrowb))
		info = 10;
	ប្រសិនបើ (lda < MAX(1, nrowa))
		info = 8;
	ប្រសិនបើ (k < 0)
		info = 5;
	ប្រសិនបើ (m < 0)
		info = 4;
	ប្រសិនបើ (transb < 0)
		info = 3;
	ប្រសិនបើ (transa < 0)
		info = 2;
	ប្រសិនបើ (uplo < 0)
		info = 1;

	ប្រសិនបើ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ត្រឡប់;
	}
#បើមិនដូច្នេះ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#បើមិនដូច្នេះ
	   void *valpha,
	   void *va, blasint LDA,
	   void *vb, blasint LDB, void *vbeta, void *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	ចំនួនគត់ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ប្រសិនបើ (order == CblasColMajor) {
		ប្រសិនបើ (Uplo == CblasUpper) uplo = 0;
		ប្រសិនបើ (Uplo == CblasLower) uplo = 1;

		ប្រសិនបើ (TransA == CblasNoTrans)
			transa = 0;
		ប្រសិនបើ (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ប្រសិនបើ (TransA == CblasConjNoTrans)
			transa = 0;
		ប្រសិនបើ (TransA == CblasConjTrans)
			transa = 1;
#បើមិនដូច្នេះ
		ប្រសិនបើ (TransA == CblasConjNoTrans)
			transa = 2;
		ប្រសិនបើ (TransA == CblasConjTrans)
			transa = 3;
#endif
		ប្រសិនបើ (TransB == CblasNoTrans)
			transb = 0;
		ប្រសិនបើ (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ប្រសិនបើ (TransB == CblasConjNoTrans)
			transb = 0;
		ប្រសិនបើ (TransB == CblasConjTrans)
			transb = 1;
#បើមិនដូច្នេះ
		ប្រសិនបើ (TransB == CblasConjNoTrans)
			transb = 2;
		ប្រសិនបើ (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ប្រសិនបើ (transa & 1) nrowa = k;
		nrowb = k;
		ប្រសិនបើ (transb & 1) nrowb = m;

		ប្រសិនបើ (ldc < MAX(1, m))
			info = 13;
		ប្រសិនបើ (ldb < MAX(1, nrowb))
			info = 10;
		ប្រសិនបើ (lda < MAX(1, nrowa))
			info = 8;
		ប្រសិនបើ (k < 0)
			info = 5;
		ប្រសិនបើ (m < 0)
			info = 4;
		ប្រសិនបើ (transb < 0)
			info = 3;
		ប្រសិនបើ (transa < 0)
			info = 2;
		ប្រសិនបើ (uplo < 0)
			info = 1;
	}

	ប្រសិនបើ (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		ប្រសិនបើ (Uplo == CblasUpper) uplo = 1;
		ប្រសិនបើ (Uplo == CblasLower) uplo = 0;

		ប្រសិនបើ (TransB == CblasNoTrans)
			transa = 0;
		ប្រសិនបើ (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ប្រសិនបើ (TransB == CblasConjNoTrans)
			transa = 0;
		ប្រសិនបើ (TransB == CblasConjTrans)
			transa = 1;
#បើមិនដូច្នេះ
		ប្រសិនបើ (TransB == CblasConjNoTrans)
			transa = 2;
		ប្រសិនបើ (TransB == CblasConjTrans)
			transa = 3;
#endif
		ប្រសិនបើ (TransA == CblasNoTrans)
			transb = 0;
		ប្រសិនបើ (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ប្រសិនបើ (TransA == CblasConjNoTrans)
			transb = 0;
		ប្រសិនបើ (TransA == CblasConjTrans)
			transb = 1;
#បើមិនដូច្នេះ
		ប្រសិនបើ (TransA == CblasConjNoTrans)
			transb = 2;
		ប្រសិនបើ (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ប្រសិនបើ (transa & 1) ncola = k;
		ncolb = k;
		ប្រសិនបើ (transb & 1) ncolb = m;

		ប្រសិនបើ (ldc < MAX(1,m))
			info = 13;
		ប្រសិនបើ (ldb < MAX(1, ncolb))
			info = 8;
		ប្រសិនបើ (lda < MAX(1, ncola))
			info = 10;
		ប្រសិនបើ (k < 0)
			info = 5;
		ប្រសិនបើ (m < 0)
			info = 4;
		ប្រសិនបើ (transb < 0)
			info = 2;
		ប្រសិនបើ (transa < 0)
			info = 3;
		ប្រសិនបើ (uplo < 0)
			info = 1;
	}

	ប្រសិនបើ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ត្រឡប់;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#បើមិនដូច្នេះ
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	dviguu_na MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	ចំនួនគត់ mode = BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
	ចំនួនគត់ mode = BLAS_DOUBLE  | BLAS_REAL;
#បើមិនដូច្នេះ
	ចំនួនគត់ mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#បើមិនដូច្នេះ
#ifdef XDOUBLE
	ចំនួនគត់ mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
	ចំនួនគត់ mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#បើមិនដូច្នេះ
	ចំនួនគត់ mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	static ចំនួនគត់ (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#បើមិនដូច្នេះ
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

	ប្រសិនបើ (m == 0)
		ត្រឡប់;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (void *)a;
	args.b = (void *)b;
	args.c = (void *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#ប្រសិនបើ defined(COMPLEX)
	args.alpha = (void *)alpha;
	args.beta  = (void *)beta;
#បើមិនដូច្នេះ
	args.alpha = (void *)&alpha;
	args.beta  = (void *)&beta;
#endif

#ifdef SMP
	args.common   = NULL;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	ប្រសិនបើ (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ត្រឡប់;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (dviguu_na)(m + 1) * (dviguu_na)m * (dviguu_na)k;
	ប្រសិនបើ (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} បើមិនដូច្នេះ {
		args.nthreads = num_cpu_avail(3);
	}

	ប្រសិនបើ (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
	} បើមិនដូច្នេះ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, NULL, NULL, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	ត្រឡប់;
}
