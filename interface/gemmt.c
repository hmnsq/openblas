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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#ନଚେତଯଦି defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#ଅନ୍ୟଥା
#define ERROR_NAME "SGEMMTR"
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#ନଚେତଯଦି defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#ଅନ୍ୟଥା
#define ERROR_NAME "SGEMMT "
#endif
#endif
#ଅନ୍ୟଥା
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMMTR"
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, ଅକ୍ଷର *TRANSA, ଅକ୍ଷର *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା transa, transb, uplo;
	blasint info;

	ଅକ୍ଷର transA, transB, Uplo;
	blasint nrowa, nrowb;
#ଯଦି !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ଯଦି defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ଅନ୍ୟଥା
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

	ଯଦି (transA == 'N')
		transa = 0;
	ଯଦି (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	ଯଦି (transA == 'R')
		transa = 0;
	ଯଦି (transA == 'C')
		transa = 1;
#ଅନ୍ୟଥା
	ଯଦି (transA == 'R')
		transa = 2;
	ଯଦି (transA == 'C')
		transa = 3;
#endif

	ଯଦି (transB == 'N')
		transb = 0;
	ଯଦି (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	ଯଦି (transB == 'R')
		transb = 0;
	ଯଦି (transB == 'C')
		transb = 1;
#ଅନ୍ୟଥା
	ଯଦି (transB == 'R')
		transb = 2;
	ଯଦି (transB == 'C')
		transb = 3;
#endif

	ଯଦି (Uplo == 'U')
		uplo = 0;
	ଯଦି (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ଯଦି (transa & 1) nrowa = k;
	nrowb = k;
	ଯଦି (transb & 1) nrowb = m;

	info = 0;

	ଯଦି (ldc < MAX(1, m))
		info = 13;
	ଯଦି (ldb < MAX(1, nrowb))
		info = 10;
	ଯଦି (lda < MAX(1, nrowa))
		info = 8;
	ଯଦି (k < 0)
		info = 5;
	ଯଦି (m < 0)
		info = 4;
	ଯଦି (transb < 0)
		info = 3;
	ଯଦି (transa < 0)
		info = 2;
	ଯଦି (uplo < 0)
		info = 1;

	ଯଦି (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ଆକାର(ERROR_NAME));
		ଫେରନ୍ତୁ;
	}
#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_UPLO Uplo,
	   ଗଣନା CBLAS_TRANSPOSE TransA, ଗଣନା CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *valpha,
	   ଶୂନ୍ୟ *va, blasint LDA,
	   ଶୂନ୍ୟ *vb, blasint LDB, ଶୂନ୍ୟ *vbeta, ଶୂନ୍ୟ *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ଯଦି (order == CblasColMajor) {
		ଯଦି (Uplo == CblasUpper) uplo = 0;
		ଯଦି (Uplo == CblasLower) uplo = 1;

		ଯଦି (TransA == CblasNoTrans)
			transa = 0;
		ଯଦି (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ଯଦି (TransA == CblasConjNoTrans)
			transa = 0;
		ଯଦି (TransA == CblasConjTrans)
			transa = 1;
#ଅନ୍ୟଥା
		ଯଦି (TransA == CblasConjNoTrans)
			transa = 2;
		ଯଦି (TransA == CblasConjTrans)
			transa = 3;
#endif
		ଯଦି (TransB == CblasNoTrans)
			transb = 0;
		ଯଦି (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ଯଦି (TransB == CblasConjNoTrans)
			transb = 0;
		ଯଦି (TransB == CblasConjTrans)
			transb = 1;
#ଅନ୍ୟଥା
		ଯଦି (TransB == CblasConjNoTrans)
			transb = 2;
		ଯଦି (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (ଶୂନ୍ୟ *)A;
		b = (ଶୂନ୍ୟ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ଯଦି (transa & 1) nrowa = k;
		nrowb = k;
		ଯଦି (transb & 1) nrowb = m;

		ଯଦି (ldc < MAX(1, m))
			info = 13;
		ଯଦି (ldb < MAX(1, nrowb))
			info = 10;
		ଯଦି (lda < MAX(1, nrowa))
			info = 8;
		ଯଦି (k < 0)
			info = 5;
		ଯଦି (m < 0)
			info = 4;
		ଯଦି (transb < 0)
			info = 3;
		ଯଦି (transa < 0)
			info = 2;
		ଯଦି (uplo < 0)
			info = 1;
	}

	ଯଦି (order == CblasRowMajor) {

		a = (ଶୂନ୍ୟ *)B;
		b = (ଶୂନ୍ୟ *)A;

		lda = LDB;
		ldb = LDA;

		ଯଦି (Uplo == CblasUpper) uplo = 1;
		ଯଦି (Uplo == CblasLower) uplo = 0;

		ଯଦି (TransB == CblasNoTrans)
			transa = 0;
		ଯଦି (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ଯଦି (TransB == CblasConjNoTrans)
			transa = 0;
		ଯଦି (TransB == CblasConjTrans)
			transa = 1;
#ଅନ୍ୟଥା
		ଯଦି (TransB == CblasConjNoTrans)
			transa = 2;
		ଯଦି (TransB == CblasConjTrans)
			transa = 3;
#endif
		ଯଦି (TransA == CblasNoTrans)
			transb = 0;
		ଯଦି (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ଯଦି (TransA == CblasConjNoTrans)
			transb = 0;
		ଯଦି (TransA == CblasConjTrans)
			transb = 1;
#ଅନ୍ୟଥା
		ଯଦି (TransA == CblasConjNoTrans)
			transb = 2;
		ଯଦି (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ଯଦି (transa & 1) ncola = k;
		ncolb = k;
		ଯଦି (transb & 1) ncolb = m;

		ଯଦି (ldc < MAX(1,m))
			info = 13;
		ଯଦି (ldb < MAX(1, ncolb))
			info = 8;
		ଯଦି (lda < MAX(1, ncola))
			info = 10;
		ଯଦି (k < 0)
			info = 5;
		ଯଦି (m < 0)
			info = 4;
		ଯଦି (transb < 0)
			info = 2;
		ଯଦି (transa < 0)
			info = 3;
		ଯଦି (uplo < 0)
			info = 1;
	}

	ଯଦି (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ଆକାର(ERROR_NAME));
		ଫେରନ୍ତୁ;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#ଅନ୍ୟଥା
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	ଦ୍ୱିଗୁଣ MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode = BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode = BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ଅନ୍ୟଥା
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

	ଯଦି (m == 0)
		ଫେରନ୍ତୁ;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (ଶୂନ୍ୟ *)a;
	args.b = (ଶୂନ୍ୟ *)b;
	args.c = (ଶୂନ୍ୟ *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#ଯଦି defined(COMPLEX)
	args.alpha = (ଶୂନ୍ୟ *)alpha;
	args.beta  = (ଶୂନ୍ୟ *)beta;
#ଅନ୍ୟଥା
	args.alpha = (ଶୂନ୍ୟ *)&alpha;
	args.beta  = (ଶୂନ୍ୟ *)&beta;
#endif

#ifdef SMP
	args.common   = ଖାଲି;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	ଯଦି (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, ଆକାର(ERROR_NAME));
		ଫେରନ୍ତୁ;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (ଦ୍ୱିଗୁଣ)(m + 1) * (ଦ୍ୱିଗୁଣ)m * (ଦ୍ୱିଗୁଣ)k;
	ଯଦି (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} ଅନ୍ୟଥା {
		args.nthreads = ଅଙ୍କ_cpu_avail(3);
	}

	ଯଦି (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP
	} ଅନ୍ୟଥା {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, ଖାଲି, ଖାଲି, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	ଫେରନ୍ତୁ;
}
