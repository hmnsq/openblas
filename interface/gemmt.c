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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <stdlib.h>
#அந்தர்பாவ "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#அல்லதெனில் defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#இல்லையெனில்
#define ERROR_NAME "SGEMMTR"
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#அல்லதெனில் defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#இல்லையெனில்
#define ERROR_NAME "SGEMMT "
#endif
#endif
#இல்லையெனில்
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#இல்லையெனில்
#define ERROR_NAME "CGEMMTR"
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#இல்லையெனில்
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, எழுத்து *TRANSA, எழுத்து *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	முழு transa, transb, uplo;
	blasint info;

	எழுத்து transA, transB, Uplo;
	blasint nrowa, nrowb;
#எனில் !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#எனில் defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#இல்லையெனில்
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

	எனில் (transA == 'N')
		transa = 0;
	எனில் (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	எனில் (transA == 'R')
		transa = 0;
	எனில் (transA == 'C')
		transa = 1;
#இல்லையெனில்
	எனில் (transA == 'R')
		transa = 2;
	எனில் (transA == 'C')
		transa = 3;
#endif

	எனில் (transB == 'N')
		transb = 0;
	எனில் (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	எனில் (transB == 'R')
		transb = 0;
	எனில் (transB == 'C')
		transb = 1;
#இல்லையெனில்
	எனில் (transB == 'R')
		transb = 2;
	எனில் (transB == 'C')
		transb = 3;
#endif

	எனில் (Uplo == 'U')
		uplo = 0;
	எனில் (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	எனில் (transa & 1) nrowa = k;
	nrowb = k;
	எனில் (transb & 1) nrowb = m;

	info = 0;

	எனில் (ldc < MAX(1, m))
		info = 13;
	எனில் (ldb < MAX(1, nrowb))
		info = 10;
	எனில் (lda < MAX(1, nrowa))
		info = 8;
	எனில் (k < 0)
		info = 5;
	எனில் (m < 0)
		info = 4;
	எனில் (transb < 0)
		info = 3;
	எனில் (transa < 0)
		info = 2;
	எனில் (uplo < 0)
		info = 1;

	எனில் (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, அளவு(ERROR_NAME));
		திரும்பு;
	}
#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_UPLO Uplo,
	   எண்ணிக்கை CBLAS_TRANSPOSE TransA, எண்ணிக்கை CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#இல்லையெனில்
	   வெறுமை *valpha,
	   வெறுமை *va, blasint LDA,
	   வெறுமை *vb, blasint LDB, வெறுமை *vbeta, வெறுமை *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	முழு transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	எனில் (order == CblasColMajor) {
		எனில் (Uplo == CblasUpper) uplo = 0;
		எனில் (Uplo == CblasLower) uplo = 1;

		எனில் (TransA == CblasNoTrans)
			transa = 0;
		எனில் (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		எனில் (TransA == CblasConjNoTrans)
			transa = 0;
		எனில் (TransA == CblasConjTrans)
			transa = 1;
#இல்லையெனில்
		எனில் (TransA == CblasConjNoTrans)
			transa = 2;
		எனில் (TransA == CblasConjTrans)
			transa = 3;
#endif
		எனில் (TransB == CblasNoTrans)
			transb = 0;
		எனில் (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		எனில் (TransB == CblasConjNoTrans)
			transb = 0;
		எனில் (TransB == CblasConjTrans)
			transb = 1;
#இல்லையெனில்
		எனில் (TransB == CblasConjNoTrans)
			transb = 2;
		எனில் (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (வெறுமை *)A;
		b = (வெறுமை *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		எனில் (transa & 1) nrowa = k;
		nrowb = k;
		எனில் (transb & 1) nrowb = m;

		எனில் (ldc < MAX(1, m))
			info = 13;
		எனில் (ldb < MAX(1, nrowb))
			info = 10;
		எனில் (lda < MAX(1, nrowa))
			info = 8;
		எனில் (k < 0)
			info = 5;
		எனில் (m < 0)
			info = 4;
		எனில் (transb < 0)
			info = 3;
		எனில் (transa < 0)
			info = 2;
		எனில் (uplo < 0)
			info = 1;
	}

	எனில் (order == CblasRowMajor) {

		a = (வெறுமை *)B;
		b = (வெறுமை *)A;

		lda = LDB;
		ldb = LDA;

		எனில் (Uplo == CblasUpper) uplo = 1;
		எனில் (Uplo == CblasLower) uplo = 0;

		எனில் (TransB == CblasNoTrans)
			transa = 0;
		எனில் (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		எனில் (TransB == CblasConjNoTrans)
			transa = 0;
		எனில் (TransB == CblasConjTrans)
			transa = 1;
#இல்லையெனில்
		எனில் (TransB == CblasConjNoTrans)
			transa = 2;
		எனில் (TransB == CblasConjTrans)
			transa = 3;
#endif
		எனில் (TransA == CblasNoTrans)
			transb = 0;
		எனில் (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		எனில் (TransA == CblasConjNoTrans)
			transb = 0;
		எனில் (TransA == CblasConjTrans)
			transb = 1;
#இல்லையெனில்
		எனில் (TransA == CblasConjNoTrans)
			transb = 2;
		எனில் (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		எனில் (transa & 1) ncola = k;
		ncolb = k;
		எனில் (transb & 1) ncolb = m;

		எனில் (ldc < MAX(1,m))
			info = 13;
		எனில் (ldb < MAX(1, ncolb))
			info = 8;
		எனில் (lda < MAX(1, ncola))
			info = 10;
		எனில் (k < 0)
			info = 5;
		எனில் (m < 0)
			info = 4;
		எனில் (transb < 0)
			info = 2;
		எனில் (transa < 0)
			info = 3;
		எனில் (uplo < 0)
			info = 1;
	}

	எனில் (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, அளவு(ERROR_NAME));
		திரும்பு;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#இல்லையெனில்
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	இரட்டை MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	முழு mode = BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
	முழு mode = BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
	முழு mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
	முழு mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
	முழு mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
	முழு mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	நிலையான முழு (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#இல்லையெனில்
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

	எனில் (m == 0)
		திரும்பு;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (வெறுமை *)a;
	args.b = (வெறுமை *)b;
	args.c = (வெறுமை *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#எனில் defined(COMPLEX)
	args.alpha = (வெறுமை *)alpha;
	args.beta  = (வெறுமை *)beta;
#இல்லையெனில்
	args.alpha = (வெறுமை *)&alpha;
	args.beta  = (வெறுமை *)&beta;
#endif

#ifdef SMP
	args.common   = வெற்று;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	எனில் (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, அளவு(ERROR_NAME));
		திரும்பு;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (இரட்டை)(m + 1) * (இரட்டை)m * (இரட்டை)k;
	எனில் (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} இல்லையெனில் {
		args.nthreads = இலக்கம்_cpu_avail(3);
	}

	எனில் (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, வெற்று, வெற்று, sa, sb, 0);

#ifdef SMP
	} இல்லையெனில் {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, வெற்று, வெற்று, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	திரும்பு;
}
