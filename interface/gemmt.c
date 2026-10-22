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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <stdlib.h>
#ಅನ್ತರ್ಭಾವ "common.h"

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMTR"
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DGEMMTR"
#ನೋಚೇದ್ಯದಿ defined(BFLOAT16)
#define ERROR_NAME "SBGEMMTR"
#ಅನ್ಯಥಾ
#define ERROR_NAME "SGEMMTR"
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "QGEMMT "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DGEMMT "
#ನೋಚೇದ್ಯದಿ defined(BFLOAT16)
#define ERROR_NAME "SBGEMMT "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SGEMMT "
#endif
#endif
#ಅನ್ಯಥಾ
#define SMP_THRESHOLD_MIN 8192.0
#ifdef RNAME
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMTR"
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMMTR"
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMMTR"
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMMT "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMMT "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMMT "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, ವರ್ಣ *TRANSA, ವರ್ಣ *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ಪೂರ್ಣಾಙ್ಕ transa, transb, uplo;
	blasint info;

	ವರ್ಣ transA, transB, Uplo;
	blasint nrowa, nrowb;
#ಯದಿ !defined(COMPLEX)
	FLOAT alpha, beta;
#endif

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ಯದಿ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ಅನ್ಯಥಾ
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

	ಯದಿ (transA == 'N')
		transa = 0;
	ಯದಿ (transA == 'T')
		transa = 1;
#ifndef COMPLEX
	ಯದಿ (transA == 'R')
		transa = 0;
	ಯದಿ (transA == 'C')
		transa = 1;
#ಅನ್ಯಥಾ
	ಯದಿ (transA == 'R')
		transa = 2;
	ಯದಿ (transA == 'C')
		transa = 3;
#endif

	ಯದಿ (transB == 'N')
		transb = 0;
	ಯದಿ (transB == 'T')
		transb = 1;
#ifndef COMPLEX
	ಯದಿ (transB == 'R')
		transb = 0;
	ಯದಿ (transB == 'C')
		transb = 1;
#ಅನ್ಯಥಾ
	ಯದಿ (transB == 'R')
		transb = 2;
	ಯದಿ (transB == 'C')
		transb = 3;
#endif

	ಯದಿ (Uplo == 'U')
		uplo = 0;
	ಯದಿ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ಯದಿ (transa & 1) nrowa = k;
	nrowb = k;
	ಯದಿ (transb & 1) nrowb = m;

	info = 0;

	ಯದಿ (ldc < MAX(1, m))
		info = 13;
	ಯದಿ (ldb < MAX(1, nrowb))
		info = 10;
	ಯದಿ (lda < MAX(1, nrowa))
		info = 8;
	ಯದಿ (k < 0)
		info = 5;
	ಯದಿ (m < 0)
		info = 4;
	ಯದಿ (transb < 0)
		info = 3;
	ಯದಿ (transa < 0)
		info = 2;
	ಯದಿ (uplo < 0)
		info = 1;

	ಯದಿ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
		ಪ್ರತ್ಯಾಗಮ;
	}
#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_UPLO Uplo,
	   ಗಣನಾ CBLAS_TRANSPOSE TransA, ಗಣನಾ CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint LDA,
	   ಶೂನ್ಯ *vb, blasint LDB, ಶೂನ್ಯ *vbeta, ಶೂನ್ಯ *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#endif
	ಪೂರ್ಣಾಙ್ಕ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ಯದಿ (order == CblasColMajor) {
		ಯದಿ (Uplo == CblasUpper) uplo = 0;
		ಯದಿ (Uplo == CblasLower) uplo = 1;

		ಯದಿ (TransA == CblasNoTrans)
			transa = 0;
		ಯದಿ (TransA == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ಯದಿ (TransA == CblasConjNoTrans)
			transa = 0;
		ಯದಿ (TransA == CblasConjTrans)
			transa = 1;
#ಅನ್ಯಥಾ
		ಯದಿ (TransA == CblasConjNoTrans)
			transa = 2;
		ಯದಿ (TransA == CblasConjTrans)
			transa = 3;
#endif
		ಯದಿ (TransB == CblasNoTrans)
			transb = 0;
		ಯದಿ (TransB == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ಯದಿ (TransB == CblasConjNoTrans)
			transb = 0;
		ಯದಿ (TransB == CblasConjTrans)
			transb = 1;
#ಅನ್ಯಥಾ
		ಯದಿ (TransB == CblasConjNoTrans)
			transb = 2;
		ಯದಿ (TransB == CblasConjTrans)
			transb = 3;
#endif

		a = (ಶೂನ್ಯ *)A;
		b = (ಶೂನ್ಯ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ಯದಿ (transa & 1) nrowa = k;
		nrowb = k;
		ಯದಿ (transb & 1) nrowb = m;

		ಯದಿ (ldc < MAX(1, m))
			info = 13;
		ಯದಿ (ldb < MAX(1, nrowb))
			info = 10;
		ಯದಿ (lda < MAX(1, nrowa))
			info = 8;
		ಯದಿ (k < 0)
			info = 5;
		ಯದಿ (m < 0)
			info = 4;
		ಯದಿ (transb < 0)
			info = 3;
		ಯದಿ (transa < 0)
			info = 2;
		ಯದಿ (uplo < 0)
			info = 1;
	}

	ಯದಿ (order == CblasRowMajor) {

		a = (ಶೂನ್ಯ *)B;
		b = (ಶೂನ್ಯ *)A;

		lda = LDB;
		ldb = LDA;

		ಯದಿ (Uplo == CblasUpper) uplo = 1;
		ಯದಿ (Uplo == CblasLower) uplo = 0;

		ಯದಿ (TransB == CblasNoTrans)
			transa = 0;
		ಯದಿ (TransB == CblasTrans)
			transa = 1;
#ifndef COMPLEX
		ಯದಿ (TransB == CblasConjNoTrans)
			transa = 0;
		ಯದಿ (TransB == CblasConjTrans)
			transa = 1;
#ಅನ್ಯಥಾ
		ಯದಿ (TransB == CblasConjNoTrans)
			transa = 2;
		ಯದಿ (TransB == CblasConjTrans)
			transa = 3;
#endif
		ಯದಿ (TransA == CblasNoTrans)
			transb = 0;
		ಯದಿ (TransA == CblasTrans)
			transb = 1;
#ifndef COMPLEX
		ಯದಿ (TransA == CblasConjNoTrans)
			transb = 0;
		ಯದಿ (TransA == CblasConjTrans)
			transb = 1;
#ಅನ್ಯಥಾ
		ಯದಿ (TransA == CblasConjNoTrans)
			transb = 2;
		ಯದಿ (TransA == CblasConjTrans)
			transb = 3;
#endif

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ಯದಿ (transa & 1) ncola = k;
		ncolb = k;
		ಯದಿ (transb & 1) ncolb = m;

		ಯದಿ (ldc < MAX(1,m))
			info = 13;
		ಯದಿ (ldb < MAX(1, ncolb))
			info = 8;
		ಯದಿ (lda < MAX(1, ncola))
			info = 10;
		ಯದಿ (k < 0)
			info = 5;
		ಯದಿ (m < 0)
			info = 4;
		ಯದಿ (transb < 0)
			info = 2;
		ಯದಿ (transa < 0)
			info = 3;
		ಯದಿ (uplo < 0)
			info = 1;
	}

	ಯದಿ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
		ಪ್ರತ್ಯಾಗಮ;
	}
#endif
#ifndef COMPLEX
#define GEMMT_TABLE_INDEX ((uplo << 2) | (transb << 1) | transa)
#ಅನ್ಯಥಾ
#define GEMMT_TABLE_INDEX ((uplo << 4) | (transb << 2) | transa)
#endif

	blas_arg_t args;
	FLOAT *buffer, *sa, *sb;

#ifdef SMP
	ದ್ವಿಗುಣ MMK;
#ifndef COMPLEX
#ifdef XDOUBLE
	ಪೂರ್ಣಾಙ್ಕ mode = BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
	ಪೂರ್ಣಾಙ್ಕ mode = BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
	ಪೂರ್ಣಾಙ್ಕ mode = BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
	ಪೂರ್ಣಾಙ್ಕ mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
	ಪೂರ್ಣಾಙ್ಕ mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
	ಪೂರ್ಣಾಙ್ಕ mode = BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ಅನ್ಯಥಾ
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

	ಯದಿ (m == 0)
		ಪ್ರತ್ಯಾಗಮ;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (ಶೂನ್ಯ *)a;
	args.b = (ಶೂನ್ಯ *)b;
	args.c = (ಶೂನ್ಯ *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#ಯದಿ defined(COMPLEX)
	args.alpha = (ಶೂನ್ಯ *)alpha;
	args.beta  = (ಶೂನ್ಯ *)beta;
#ಅನ್ಯಥಾ
	args.alpha = (ಶೂನ್ಯ *)&alpha;
	args.beta  = (ಶೂನ್ಯ *)&beta;
#endif

#ifdef SMP
	args.common   = ರಿಕ್ತ;
	args.nthreads = 1;
#endif

	buffer = (FLOAT *)blas_memory_alloc(0);
	ಯದಿ (!buffer) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
		ಪ್ರತ್ಯಾಗಮ;
	}

	sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
	MMK = (ದ್ವಿಗುಣ)(m + 1) * (ದ್ವಿಗುಣ)m * (ದ್ವಿಗುಣ)k;
	ಯದಿ (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} ಅನ್ಯಥಾ {
		args.nthreads = ಅಂಕೆ_cpu_avail(3);
	}

	ಯದಿ (args.nthreads == 1) {
#endif

		(gemmt[GEMMT_TABLE_INDEX])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ifdef SMP
	} ಅನ್ಯಥಾ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, ರಿಕ್ತ, ರಿಕ್ತ, gemmt[GEMMT_TABLE_INDEX], sa, sb, args.nthreads);
	}
#endif

	blas_memory_free(buffer);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	ಪ್ರತ್ಯಾಗಮ;
}
