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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <stdlib.h>
#અન્તર્ભાવ "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#અથવાજો defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#નહીંતર
#define ERROR_NAME "SGEMMTR"
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#અથવાજો defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#નહીંતર
#define ERROR_NAME "SGEMMT "
#endif
#endif
#નહીંતર
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#નહીંતર
#define ERROR_NAME "CGEMMTR"
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#નહીંતર
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, અક્ષર *TRANSA, અક્ષર *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	પૂર્ણાંક transa, transb, uplo;
	blasint info;

	અક્ષર transA, transB, Uplo;
	blasint nrowa, nrowb;
#જો !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#જો defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#નહીંતર
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

	જો (transA == 'N')
		transa = 0;
	જો (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	જો (transA == 'R')
		transa = 0;
	જો (transA == 'C')
		transa = 1;
#નહીંતર
	જો (transA == 'R')
		transa = 2;
	જો (transA == 'C')
		transa = 3;
#endif

	જો (transB == 'N')
		transb = 0;
	જો (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	જો (transB == 'R')
		transb = 0;
	જો (transB == 'C')
		transb = 1;
#નહીંતર
	જો (transB == 'R')
		transb = 2;
	જો (transB == 'C')
		transb = 3;
#endif

	જો (Uplo == 'U')
		uplo = 0;
	જો (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	જો (transa & 1) nrowa = k;
	nrowb = k;
	જો (transb & 1) nrowb = m;

	info = 0;

	જો (ldc < MAX(1, m))
		info = 13;
	જો (ldb < MAX(1, nrowb))
		info = 10;
	જો (lda < MAX(1, nrowa))
		info = 8;
	જો (k < 0)
		info = 5;
	જો (m < 0)
		info = 4;
	જો (transb < 0)
		info = 3;
	જો (transa < 0)
		info = 2;
	જો (uplo < 0)
		info = 1;

	જો (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, કદ(ERROR_NAME));
		પાછા;
	}
#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_UPLO Uplo,
	   ગણના CBLAS_TRANSPOSE TransA, ગણના CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#નહીંતર
	   શૂન્ય *valpha,
	   શૂન્ય *va, blasint LDA,
	   શૂન્ય *vb, blasint LDB, શૂન્ય *vbeta, શૂન્ય *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	પૂર્ણાંક transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	જો (order == CblasColMajor) {
		જો (Uplo == CblasUpper) uplo = 0;
		જો (Uplo == CblasLower) uplo = 1;

		જો (TransA == CblasNoTrans)
			transa = 0;
		જો (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		જો (TransA == CblasConjNoTrans)
			transa = 0;
		જો (TransA == CblasConjTrans)
			transa = 1;
#નહીંતર
		જો (TransA == CblasConjNoTrans)
			transa = 2;
		જો (TransA == CblasConjTrans)
			transa = 3;
#endif
		જો (TransB == CblasNoTrans)
			transb = 0;
		જો (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		જો (TransB == CblasConjNoTrans)
			transb = 0;
		જો (TransB == CblasConjTrans)
			transb = 1;
#નહીંતર
		જો (TransB == CblasConjNoTrans)
			transb = 2;
		જો (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (શૂન્ય *)A;
		b = (શૂન્ય *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		જો (transa & 1) nrowa = k;
		nrowb = k;
		જો (transb & 1) nrowb = m;

		જો (ldc < MAX(1, m))
			info = 13;
		જો (ldb < MAX(1, nrowb))
			info = 10;
		જો (lda < MAX(1, nrowa))
			info = 8;
		જો (k < 0)
			info = 5;
		જો (m < 0)
			info = 4;
		જો (transb < 0)
			info = 3;
		જો (transa < 0)
			info = 2;
		જો (uplo < 0)
			info = 1;
	}

	જો (order == CblasRowMajor) {

		a = (શૂન્ય *)B;
		b = (શૂન્ય *)A;

		lda = LDB;
		ldb = LDA;

		જો (Uplo == CblasUpper) uplo = 1;
		જો (Uplo == CblasLower) uplo = 0;

		જો (TransB == CblasNoTrans)
			transa = 0;
		જો (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		જો (TransB == CblasConjNoTrans)
			transa = 0;
		જો (TransB == CblasConjTrans)
			transa = 1;
#નહીંતર
		જો (TransB == CblasConjNoTrans)
			transa = 2;
		જો (TransB == CblasConjTrans)
			transa = 3;
#endif
		જો (TransA == CblasNoTrans)
			transb = 0;
		જો (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		જો (TransA == CblasConjNoTrans)
			transb = 0;
		જો (TransA == CblasConjTrans)
			transb = 1;
#નહીંતર
		જો (TransA == CblasConjNoTrans)
			transb = 2;
		જો (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		જો (transa & 1) ncola = k;
		ncolb = k;
		જો (transb & 1) ncolb = m;

		જો (ldc < MAX(1,m))
			info = 13;
		જો (ldb < MAX(1, ncolb))
			info = 8;
		જો (lda < MAX(1, ncola))
			info = 10;
		જો (k < 0)
			info = 5;
		જો (m < 0)
			info = 4;
		જો (transb < 0)
			info = 2;
		જો (transa < 0)
			info = 3;
		જો (uplo < 0)
			info = 1;
	}

	જો (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, કદ(ERROR_NAME));
		પાછા;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#નહીંતર
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	દ્વિગુણ MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	પૂર્ણાંક mode = BLAS_XDOUBLE | BLAS_REAL;
#અથવાજો defined(DOUBLE)
	પૂર્ણાંક mode = BLAS_DOUBLE  | BLAS_REAL;
#નહીંતર
	પૂર્ણાંક mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#નહીંતર
#ifdef XDOUBLE
	પૂર્ણાંક mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
	પૂર્ણાંક mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
	પૂર્ણાંક mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	સ્થિર પૂર્ણાંક (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#નહીંતર
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

	જો (m == 0)
		પાછા;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (શૂન્ય *)a;
	args.b = (શૂન્ય *)b;
	args.c = (શૂન્ય *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#જો defined(COMPLEX)
	args.alpha = (શૂન્ય *)alpha;
	args.beta  = (શૂન્ય *)beta;
#નહીંતર
	args.alpha = (શૂન્ય *)&alpha;
	args.beta  = (શૂન્ય *)&beta;
#endif

#ifdef SMP
	args.common   = ખાલી;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	જો (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, કદ(ERROR_NAME));
		પાછા;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (દ્વિગુણ)(m + 1) * (દ્વિગુણ)m * (દ્વિગુણ)k;
	જો (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} નહીંતર {
		args.nthreads = અંક_cpu_avail(3);
	}

	જો (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, ખાલી, ખાલી, sa, sb, 0);

#ifdef SMP
	} નહીંતર {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, ખાલી, ખાલી, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	પાછા;
}
