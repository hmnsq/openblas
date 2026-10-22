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
#ᬮᬶᬬ
#define ERROR_NAME "SGEMMTR"
#endif
#ᬮᬶᬬ
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#elif defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#ᬮᬶᬬ
#define ERROR_NAME "SGEMMT "
#endif
#endif
#ᬮᬶᬬ
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#ᬮᬶᬬ
#define ERROR_NAME "CGEMMTR"
#endif
#ᬮᬶᬬ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#ᬮᬶᬬ
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
	puur_naa_nka transa, transb, uplo;
	blasint info;

	var_na transA, transB, Uplo;
	blasint nrowa, nrowb;
#ᬬᬾᬦ᭄ !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ᬬᬾᬦ᭄ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ᬮᬶᬬ
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

	ᬬᬾᬦ᭄ (transA == 'N')
		transa = 0;
	ᬬᬾᬦ᭄ (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	ᬬᬾᬦ᭄ (transA == 'R')
		transa = 0;
	ᬬᬾᬦ᭄ (transA == 'C')
		transa = 1;
#ᬮᬶᬬ
	ᬬᬾᬦ᭄ (transA == 'R')
		transa = 2;
	ᬬᬾᬦ᭄ (transA == 'C')
		transa = 3;
#endif

	ᬬᬾᬦ᭄ (transB == 'N')
		transb = 0;
	ᬬᬾᬦ᭄ (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	ᬬᬾᬦ᭄ (transB == 'R')
		transb = 0;
	ᬬᬾᬦ᭄ (transB == 'C')
		transb = 1;
#ᬮᬶᬬ
	ᬬᬾᬦ᭄ (transB == 'R')
		transb = 2;
	ᬬᬾᬦ᭄ (transB == 'C')
		transb = 3;
#endif

	ᬬᬾᬦ᭄ (Uplo == 'U')
		uplo = 0;
	ᬬᬾᬦ᭄ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ᬬᬾᬦ᭄ (transa & 1) nrowa = k;
	nrowb = k;
	ᬬᬾᬦ᭄ (transb & 1) nrowb = m;

	info = 0;

	ᬬᬾᬦ᭄ (ldc < MAX(1, m))
		info = 13;
	ᬬᬾᬦ᭄ (ldb < MAX(1, nrowb))
		info = 10;
	ᬬᬾᬦ᭄ (lda < MAX(1, nrowa))
		info = 8;
	ᬬᬾᬦ᭄ (k < 0)
		info = 5;
	ᬬᬾᬦ᭄ (m < 0)
		info = 4;
	ᬬᬾᬦ᭄ (transb < 0)
		info = 3;
	ᬬᬾᬦ᭄ (transa < 0)
		info = 2;
	ᬬᬾᬦ᭄ (uplo < 0)
		info = 1;

	ᬬᬾᬦ᭄ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ᬯᬮᬶ;
	}
#ᬮᬶᬬ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ᬮᬶᬬ
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
	puur_naa_nka transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ᬬᬾᬦ᭄ (order == CblasColMajor) {
		ᬬᬾᬦ᭄ (Uplo == CblasUpper) uplo = 0;
		ᬬᬾᬦ᭄ (Uplo == CblasLower) uplo = 1;

		ᬬᬾᬦ᭄ (TransA == CblasNoTrans)
			transa = 0;
		ᬬᬾᬦ᭄ (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ᬬᬾᬦ᭄ (TransA == CblasConjNoTrans)
			transa = 0;
		ᬬᬾᬦ᭄ (TransA == CblasConjTrans)
			transa = 1;
#ᬮᬶᬬ
		ᬬᬾᬦ᭄ (TransA == CblasConjNoTrans)
			transa = 2;
		ᬬᬾᬦ᭄ (TransA == CblasConjTrans)
			transa = 3;
#endif
		ᬬᬾᬦ᭄ (TransB == CblasNoTrans)
			transb = 0;
		ᬬᬾᬦ᭄ (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ᬬᬾᬦ᭄ (TransB == CblasConjNoTrans)
			transb = 0;
		ᬬᬾᬦ᭄ (TransB == CblasConjTrans)
			transb = 1;
#ᬮᬶᬬ
		ᬬᬾᬦ᭄ (TransB == CblasConjNoTrans)
			transb = 2;
		ᬬᬾᬦ᭄ (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ᬬᬾᬦ᭄ (transa & 1) nrowa = k;
		nrowb = k;
		ᬬᬾᬦ᭄ (transb & 1) nrowb = m;

		ᬬᬾᬦ᭄ (ldc < MAX(1, m))
			info = 13;
		ᬬᬾᬦ᭄ (ldb < MAX(1, nrowb))
			info = 10;
		ᬬᬾᬦ᭄ (lda < MAX(1, nrowa))
			info = 8;
		ᬬᬾᬦ᭄ (k < 0)
			info = 5;
		ᬬᬾᬦ᭄ (m < 0)
			info = 4;
		ᬬᬾᬦ᭄ (transb < 0)
			info = 3;
		ᬬᬾᬦ᭄ (transa < 0)
			info = 2;
		ᬬᬾᬦ᭄ (uplo < 0)
			info = 1;
	}

	ᬬᬾᬦ᭄ (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		ᬬᬾᬦ᭄ (Uplo == CblasUpper) uplo = 1;
		ᬬᬾᬦ᭄ (Uplo == CblasLower) uplo = 0;

		ᬬᬾᬦ᭄ (TransB == CblasNoTrans)
			transa = 0;
		ᬬᬾᬦ᭄ (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ᬬᬾᬦ᭄ (TransB == CblasConjNoTrans)
			transa = 0;
		ᬬᬾᬦ᭄ (TransB == CblasConjTrans)
			transa = 1;
#ᬮᬶᬬ
		ᬬᬾᬦ᭄ (TransB == CblasConjNoTrans)
			transa = 2;
		ᬬᬾᬦ᭄ (TransB == CblasConjTrans)
			transa = 3;
#endif
		ᬬᬾᬦ᭄ (TransA == CblasNoTrans)
			transb = 0;
		ᬬᬾᬦ᭄ (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ᬬᬾᬦ᭄ (TransA == CblasConjNoTrans)
			transb = 0;
		ᬬᬾᬦ᭄ (TransA == CblasConjTrans)
			transb = 1;
#ᬮᬶᬬ
		ᬬᬾᬦ᭄ (TransA == CblasConjNoTrans)
			transb = 2;
		ᬬᬾᬦ᭄ (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ᬬᬾᬦ᭄ (transa & 1) ncola = k;
		ncolb = k;
		ᬬᬾᬦ᭄ (transb & 1) ncolb = m;

		ᬬᬾᬦ᭄ (ldc < MAX(1,m))
			info = 13;
		ᬬᬾᬦ᭄ (ldb < MAX(1, ncolb))
			info = 8;
		ᬬᬾᬦ᭄ (lda < MAX(1, ncola))
			info = 10;
		ᬬᬾᬦ᭄ (k < 0)
			info = 5;
		ᬬᬾᬦ᭄ (m < 0)
			info = 4;
		ᬬᬾᬦ᭄ (transb < 0)
			info = 2;
		ᬬᬾᬦ᭄ (transa < 0)
			info = 3;
		ᬬᬾᬦ᭄ (uplo < 0)
			info = 1;
	}

	ᬬᬾᬦ᭄ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ᬯᬮᬶ;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#ᬮᬶᬬ
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	dviguu_na MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	puur_naa_nka mode = BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
	puur_naa_nka mode = BLAS_DOUBLE  | BLAS_REAL;
#ᬮᬶᬬ
	puur_naa_nka mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#ᬮᬶᬬ
#ifdef XDOUBLE
	puur_naa_nka mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
	puur_naa_nka mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#ᬮᬶᬬ
	puur_naa_nka mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	static puur_naa_nka (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ᬮᬶᬬ
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

	ᬬᬾᬦ᭄ (m == 0)
		ᬯᬮᬶ;

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

#ᬬᬾᬦ᭄ defined(COMPLEX)
	args.alpha = (void *)alpha;
	args.beta  = (void *)beta;
#ᬮᬶᬬ
	args.alpha = (void *)&alpha;
	args.beta  = (void *)&beta;
#endif

#ifdef SMP
	args.common   = NULL;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	ᬬᬾᬦ᭄ (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ᬯᬮᬶ;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (dviguu_na)(m + 1) * (dviguu_na)m * (dviguu_na)k;
	ᬬᬾᬦ᭄ (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} ᬮᬶᬬ {
		args.nthreads = num_cpu_avail(3);
	}

	ᬬᬾᬦ᭄ (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
	} ᬮᬶᬬ {

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

	ᬯᬮᬶ;
}
