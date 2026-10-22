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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <stdlib.h>
#അന്തര്ഭാവ "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#അല്ലെങ്കിലെങ്കിൽ defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGEMMTR"
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#അല്ലെങ്കിലെങ്കിൽ defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGEMMT "
#endif
#endif
#അല്ലെങ്കിൽ
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#അല്ലെങ്കിൽ
#define ERROR_NAME "CGEMMTR"
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, അക്ഷരം *TRANSA, അക്ഷരം *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	പൂർണ്ണസംഖ്യ transa, transb, uplo;
	blasint info;

	അക്ഷരം transA, transB, Uplo;
	blasint nrowa, nrowb;
#എങ്കിൽ !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#എങ്കിൽ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#അല്ലെങ്കിൽ
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

	എങ്കിൽ (transA == 'N')
		transa = 0;
	എങ്കിൽ (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	എങ്കിൽ (transA == 'R')
		transa = 0;
	എങ്കിൽ (transA == 'C')
		transa = 1;
#അല്ലെങ്കിൽ
	എങ്കിൽ (transA == 'R')
		transa = 2;
	എങ്കിൽ (transA == 'C')
		transa = 3;
#endif

	എങ്കിൽ (transB == 'N')
		transb = 0;
	എങ്കിൽ (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	എങ്കിൽ (transB == 'R')
		transb = 0;
	എങ്കിൽ (transB == 'C')
		transb = 1;
#അല്ലെങ്കിൽ
	എങ്കിൽ (transB == 'R')
		transb = 2;
	എങ്കിൽ (transB == 'C')
		transb = 3;
#endif

	എങ്കിൽ (Uplo == 'U')
		uplo = 0;
	എങ്കിൽ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	എങ്കിൽ (transa & 1) nrowa = k;
	nrowb = k;
	എങ്കിൽ (transb & 1) nrowb = m;

	info = 0;

	എങ്കിൽ (ldc < MAX(1, m))
		info = 13;
	എങ്കിൽ (ldb < MAX(1, nrowb))
		info = 10;
	എങ്കിൽ (lda < MAX(1, nrowa))
		info = 8;
	എങ്കിൽ (k < 0)
		info = 5;
	എങ്കിൽ (m < 0)
		info = 4;
	എങ്കിൽ (transb < 0)
		info = 3;
	എങ്കിൽ (transa < 0)
		info = 2;
	എങ്കിൽ (uplo < 0)
		info = 1;

	എങ്കിൽ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
		തിരികെ;
	}
#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA, എണ്ണൽ CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#അല്ലെങ്കിൽ
	   ശൂന്യം *valpha,
	   ശൂന്യം *va, blasint LDA,
	   ശൂന്യം *vb, blasint LDB, ശൂന്യം *vbeta, ശൂന്യം *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	പൂർണ്ണസംഖ്യ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	എങ്കിൽ (order == CblasColMajor) {
		എങ്കിൽ (Uplo == CblasUpper) uplo = 0;
		എങ്കിൽ (Uplo == CblasLower) uplo = 1;

		എങ്കിൽ (TransA == CblasNoTrans)
			transa = 0;
		എങ്കിൽ (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		എങ്കിൽ (TransA == CblasConjNoTrans)
			transa = 0;
		എങ്കിൽ (TransA == CblasConjTrans)
			transa = 1;
#അല്ലെങ്കിൽ
		എങ്കിൽ (TransA == CblasConjNoTrans)
			transa = 2;
		എങ്കിൽ (TransA == CblasConjTrans)
			transa = 3;
#endif
		എങ്കിൽ (TransB == CblasNoTrans)
			transb = 0;
		എങ്കിൽ (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		എങ്കിൽ (TransB == CblasConjNoTrans)
			transb = 0;
		എങ്കിൽ (TransB == CblasConjTrans)
			transb = 1;
#അല്ലെങ്കിൽ
		എങ്കിൽ (TransB == CblasConjNoTrans)
			transb = 2;
		എങ്കിൽ (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (ശൂന്യം *)A;
		b = (ശൂന്യം *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		എങ്കിൽ (transa & 1) nrowa = k;
		nrowb = k;
		എങ്കിൽ (transb & 1) nrowb = m;

		എങ്കിൽ (ldc < MAX(1, m))
			info = 13;
		എങ്കിൽ (ldb < MAX(1, nrowb))
			info = 10;
		എങ്കിൽ (lda < MAX(1, nrowa))
			info = 8;
		എങ്കിൽ (k < 0)
			info = 5;
		എങ്കിൽ (m < 0)
			info = 4;
		എങ്കിൽ (transb < 0)
			info = 3;
		എങ്കിൽ (transa < 0)
			info = 2;
		എങ്കിൽ (uplo < 0)
			info = 1;
	}

	എങ്കിൽ (order == CblasRowMajor) {

		a = (ശൂന്യം *)B;
		b = (ശൂന്യം *)A;

		lda = LDB;
		ldb = LDA;

		എങ്കിൽ (Uplo == CblasUpper) uplo = 1;
		എങ്കിൽ (Uplo == CblasLower) uplo = 0;

		എങ്കിൽ (TransB == CblasNoTrans)
			transa = 0;
		എങ്കിൽ (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		എങ്കിൽ (TransB == CblasConjNoTrans)
			transa = 0;
		എങ്കിൽ (TransB == CblasConjTrans)
			transa = 1;
#അല്ലെങ്കിൽ
		എങ്കിൽ (TransB == CblasConjNoTrans)
			transa = 2;
		എങ്കിൽ (TransB == CblasConjTrans)
			transa = 3;
#endif
		എങ്കിൽ (TransA == CblasNoTrans)
			transb = 0;
		എങ്കിൽ (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		എങ്കിൽ (TransA == CblasConjNoTrans)
			transb = 0;
		എങ്കിൽ (TransA == CblasConjTrans)
			transb = 1;
#അല്ലെങ്കിൽ
		എങ്കിൽ (TransA == CblasConjNoTrans)
			transb = 2;
		എങ്കിൽ (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		എങ്കിൽ (transa & 1) ncola = k;
		ncolb = k;
		എങ്കിൽ (transb & 1) ncolb = m;

		എങ്കിൽ (ldc < MAX(1,m))
			info = 13;
		എങ്കിൽ (ldb < MAX(1, ncolb))
			info = 8;
		എങ്കിൽ (lda < MAX(1, ncola))
			info = 10;
		എങ്കിൽ (k < 0)
			info = 5;
		എങ്കിൽ (m < 0)
			info = 4;
		എങ്കിൽ (transb < 0)
			info = 2;
		എങ്കിൽ (transa < 0)
			info = 3;
		എങ്കിൽ (uplo < 0)
			info = 1;
	}

	എങ്കിൽ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
		തിരികെ;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#അല്ലെങ്കിൽ
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	ഇരട്ട MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	പൂർണ്ണസംഖ്യ mode = BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
	പൂർണ്ണസംഖ്യ mode = BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
	പൂർണ്ണസംഖ്യ mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
	പൂർണ്ണസംഖ്യ mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
	പൂർണ്ണസംഖ്യ mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
	പൂർണ്ണസംഖ്യ mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#അല്ലെങ്കിൽ
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

	എങ്കിൽ (m == 0)
		തിരികെ;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (ശൂന്യം *)a;
	args.b = (ശൂന്യം *)b;
	args.c = (ശൂന്യം *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#എങ്കിൽ defined(COMPLEX)
	args.alpha = (ശൂന്യം *)alpha;
	args.beta  = (ശൂന്യം *)beta;
#അല്ലെങ്കിൽ
	args.alpha = (ശൂന്യം *)&alpha;
	args.beta  = (ശൂന്യം *)&beta;
#endif

#ifdef SMP
	args.common   = ശൂന്യ;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	എങ്കിൽ (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
		തിരികെ;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (ഇരട്ട)(m + 1) * (ഇരട്ട)m * (ഇരട്ട)k;
	എങ്കിൽ (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} അല്ലെങ്കിൽ {
		args.nthreads = അക്കം_cpu_avail(3);
	}

	എങ്കിൽ (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#ifdef SMP
	} അല്ലെങ്കിൽ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, ശൂന്യ, ശൂന്യ, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	തിരികെ;
}
