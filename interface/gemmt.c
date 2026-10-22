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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdlib.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 65536.0
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RNAME
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGEMMTR"
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGEMMTR"
#ਜਾਂਜੇ defined(BFLOAT16)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SBGEMMTR"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGEMMTR"
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGEMMT "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGEMMT "
#ਜਾਂਜੇ defined(BFLOAT16)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SBGEMMT "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGEMMT "
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 8192.0
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RNAME
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGEMMTR"
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGEMMTR"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGEMMTR"
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGEMMT "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGEMMT "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGEMMT "
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ GEMM_MULTITHREAD_ਹੱਦ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_MULTITHREAD_ਹੱਦ 4
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, ਅੱਖਰ *TRANSA, ਅੱਖਰ *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ਪੂਰਨ transa, transb, uplo;
	blasint info;

	ਅੱਖਰ transA, transB, Uplo;
	blasint nrowa, nrowb;
#ਜੇ !defined(COMPLEX)
	FLOAT alpha, beta;
#ਜੇ_ਅੰਤ

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ਜੇ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ਨਹੀਂਤਾਂ
	alpha = *Alpha;
	beta = *Beta;
#ਜੇ_ਅੰਤ

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

	ਜੇ (transA == 'N')
		transa = 0;
	ਜੇ (transA == 'T')
		transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
	ਜੇ (transA == 'R')
		transa = 0;
	ਜੇ (transA == 'C')
		transa = 1;
#ਨਹੀਂਤਾਂ
	ਜੇ (transA == 'R')
		transa = 2;
	ਜੇ (transA == 'C')
		transa = 3;
#ਜੇ_ਅੰਤ

	ਜੇ (transB == 'N')
		transb = 0;
	ਜੇ (transB == 'T')
		transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
	ਜੇ (transB == 'R')
		transb = 0;
	ਜੇ (transB == 'C')
		transb = 1;
#ਨਹੀਂਤਾਂ
	ਜੇ (transB == 'R')
		transb = 2;
	ਜੇ (transB == 'C')
		transb = 3;
#ਜੇ_ਅੰਤ

	ਜੇ (Uplo == 'U')
		uplo = 0;
	ਜੇ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ਜੇ (transa & 1) nrowa = k;
	nrowb = k;
	ਜੇ (transb & 1) nrowb = m;

	info = 0;

	ਜੇ (ldc < MAX(1, m))
		info = 13;
	ਜੇ (ldb < MAX(1, nrowb))
		info = 10;
	ਜੇ (lda < MAX(1, nrowa))
		info = 8;
	ਜੇ (k < 0)
		info = 5;
	ਜੇ (m < 0)
		info = 4;
	ਜੇ (transb < 0)
		info = 3;
	ਜੇ (transa < 0)
		info = 2;
	ਜੇ (uplo < 0)
		info = 1;

	ਜੇ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
		ਵਾਪਸ;
	}
#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_UPLO Uplo,
	   ਗਿਣਤੀ CBLAS_TRANSPOSE TransA, ਗਿਣਤੀ CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *valpha,
	   ਖਾਲੀ *va, blasint LDA,
	   ਖਾਲੀ *vb, blasint LDB, ਖਾਲੀ *vbeta, ਖਾਲੀ *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#ਜੇ_ਅੰਤ
	ਪੂਰਨ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
		ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo = 0;
		ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo = 1;

		ਜੇ (TransA == CblasNoTrans)
			transa = 0;
		ਜੇ (TransA == CblasTrans)
			transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
		ਜੇ (TransA == CblasConjNoTrans)
			transa = 0;
		ਜੇ (TransA == CblasConjTrans)
			transa = 1;
#ਨਹੀਂਤਾਂ
		ਜੇ (TransA == CblasConjNoTrans)
			transa = 2;
		ਜੇ (TransA == CblasConjTrans)
			transa = 3;
#ਜੇ_ਅੰਤ
		ਜੇ (TransB == CblasNoTrans)
			transb = 0;
		ਜੇ (TransB == CblasTrans)
			transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
		ਜੇ (TransB == CblasConjNoTrans)
			transb = 0;
		ਜੇ (TransB == CblasConjTrans)
			transb = 1;
#ਨਹੀਂਤਾਂ
		ਜੇ (TransB == CblasConjNoTrans)
			transb = 2;
		ਜੇ (TransB == CblasConjTrans)
			transb = 3;
#ਜੇ_ਅੰਤ

		a = (ਖਾਲੀ *)A;
		b = (ਖਾਲੀ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ਜੇ (transa & 1) nrowa = k;
		nrowb = k;
		ਜੇ (transb & 1) nrowb = m;

		ਜੇ (ldc < MAX(1, m))
			info = 13;
		ਜੇ (ldb < MAX(1, nrowb))
			info = 10;
		ਜੇ (lda < MAX(1, nrowa))
			info = 8;
		ਜੇ (k < 0)
			info = 5;
		ਜੇ (m < 0)
			info = 4;
		ਜੇ (transb < 0)
			info = 3;
		ਜੇ (transa < 0)
			info = 2;
		ਜੇ (uplo < 0)
			info = 1;
	}

	ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {

		a = (ਖਾਲੀ *)B;
		b = (ਖਾਲੀ *)A;

		lda = LDB;
		ldb = LDA;

		ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo = 1;
		ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo = 0;

		ਜੇ (TransB == CblasNoTrans)
			transa = 0;
		ਜੇ (TransB == CblasTrans)
			transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
		ਜੇ (TransB == CblasConjNoTrans)
			transa = 0;
		ਜੇ (TransB == CblasConjTrans)
			transa = 1;
#ਨਹੀਂਤਾਂ
		ਜੇ (TransB == CblasConjNoTrans)
			transa = 2;
		ਜੇ (TransB == CblasConjTrans)
			transa = 3;
#ਜੇ_ਅੰਤ
		ਜੇ (TransA == CblasNoTrans)
			transb = 0;
		ਜੇ (TransA == CblasTrans)
			transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
		ਜੇ (TransA == CblasConjNoTrans)
			transb = 0;
		ਜੇ (TransA == CblasConjTrans)
			transb = 1;
#ਨਹੀਂਤਾਂ
		ਜੇ (TransA == CblasConjNoTrans)
			transb = 2;
		ਜੇ (TransA == CblasConjTrans)
			transb = 3;
#ਜੇ_ਅੰਤ

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ਜੇ (transa & 1) ncola = k;
		ncolb = k;
		ਜੇ (transb & 1) ncolb = m;

		ਜੇ (ldc < MAX(1,m))
			info = 13;
		ਜੇ (ldb < MAX(1, ncolb))
			info = 8;
		ਜੇ (lda < MAX(1, ncola))
			info = 10;
		ਜੇ (k < 0)
			info = 5;
		ਜੇ (m < 0)
			info = 4;
		ਜੇ (transb < 0)
			info = 2;
		ਜੇ (transa < 0)
			info = 3;
		ਜੇ (uplo < 0)
			info = 1;
	}

	ਜੇ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
		ਵਾਪਸ;
	}
#ਜੇ_ਅੰਤ
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMMT_ਸਾਰਣੀ_ਸੂਚਕਾਂਕ ((uplo << 2) | (transb << 1) | transa)
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMMT_ਸਾਰਣੀ_ਸੂਚਕਾਂਕ ((uplo << 4) | (transb << 2) | transa)
#ਜੇ_ਅੰਤ

	blas_arg_t ਦਲੀਲਾਂ;
	FLOAT *ਬਫਰ, *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
	ਦੁੱਗਣਾ MMK;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
	ਪੂਰਨ ਬਹੁਲਕ = BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
	ਪੂਰਨ ਬਹੁਲਕ = BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
	ਪੂਰਨ ਬਹੁਲਕ = BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
	ਪੂਰਨ ਬਹੁਲਕ = BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
	ਪੂਰਨ ਬਹੁਲਕ = BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
	ਪੂਰਨ ਬਹੁਲਕ = BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ਸਥਿਰ ਪੂਰਨ (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ਨਹੀਂਤਾਂ
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#ਜੇ_ਅੰਤ
	};

	ਜੇ (m == 0)
		ਵਾਪਸ;

	IDEBUG_ਸ਼ੁਰੂ;

	FUNCTION_PROFILE_ਸ਼ੁਰੂ();

	ਦਲੀਲਾਂ.m = m;
	ਦਲੀਲਾਂ.n = m;
	ਦਲੀਲਾਂ.k = k;

	ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
	ਦਲੀਲਾਂ.b = (ਖਾਲੀ *)b;
	ਦਲੀਲਾਂ.c = (ਖਾਲੀ *)c;

	ਦਲੀਲਾਂ.lda = lda;
	ਦਲੀਲਾਂ.ldb = ldb;
	ਦਲੀਲਾਂ.ldc = ldc;

#ਜੇ defined(COMPLEX)
	ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)alpha;
	ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)beta;
#ਨਹੀਂਤਾਂ
	ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)&alpha;
	ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)&beta;
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
	ਦਲੀਲਾਂ.common   = ਨਿਰਾ;
	ਦਲੀਲਾਂ.nthreads = 1;
#ਜੇ_ਅੰਤ

	ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(0);
	ਜੇ (!ਬਫਰ) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
		ਵਾਪਸ;
	}

	sa = (FLOAT *)((BLASLONG)ਬਫਰ + GEMM_ਵਿੱਥ_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ਮਾਪ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ਵਿੱਥ_B);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
	MMK = (ਦੁੱਗਣਾ)(m + 1) * (ਦੁੱਗਣਾ)m * (ਦੁੱਗਣਾ)k;
	ਜੇ (MMK <= (SMP_ਹੱਦ_MIN * GEMM_MULTITHREAD_ਹੱਦ)) {
		ਦਲੀਲਾਂ.nthreads = 1;
	} ਨਹੀਂਤਾਂ {
		ਦਲੀਲਾਂ.nthreads = ਅੰਕ_cpu_avail(3);
	}

	ਜੇ (ਦਲੀਲਾਂ.nthreads == 1) {
#ਜੇ_ਅੰਤ

		(gemmt[GEMMT_ਸਾਰਣੀ_ਸੂਚਕਾਂਕ])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
	} ਨਹੀਂਤਾਂ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		ਬਹੁਲਕ |= (uplo << BLAS_UPLO_ਖਿਸਕਾਓ);

		syrk_thread(ਬਹੁਲਕ, &ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, gemmt[GEMMT_ਸਾਰਣੀ_ਸੂਚਕਾਂਕ], sa, sb, ਦਲੀਲਾਂ.nthreads);
	}
#ਜੇ_ਅੰਤ

	blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

	FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE * COMPSIZE, ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.k + ਦਲੀਲਾਂ.k * ਦਲੀਲਾਂ.n + ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n / 2,
			     ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.k);

	IDEBUG_ਸਿਰਾ;

	ਵਾਪਸ;
}
