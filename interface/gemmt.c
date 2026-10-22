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

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 65536.0
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RNAME
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGEMMTR"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGEMMTR"
#ಅಥವಾಆದರೆ defined(BFLOAT16)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SBGEMMTR"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGEMMTR"
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGEMMT "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGEMMT "
#ಅಥವಾಆದರೆ defined(BFLOAT16)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SBGEMMT "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGEMMT "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 8192.0
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RNAME
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGEMMTR"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEMMTR"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEMMTR"
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGEMMT "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEMMT "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEMMT "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ GEMM_MULTITHREAD_ಮಿತಿ
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_MULTITHREAD_ಮಿತಿ 4
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *UPLO, ಅಕ್ಷರ *TRANSA, ಅಕ್ಷರ *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ಪೂರ್ಣಾಂಕ transa, transb, uplo;
	blasint info;

	ಅಕ್ಷರ transA, transB, Uplo;
	blasint nrowa, nrowb;
#ಆದರೆ !defined(COMPLEX)
	FLOAT alpha, beta;
#ಆದರೆ_ಅಂತ್ಯ

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ಆದರೆ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ಇಲ್ಲದಿದ್ದರೆ
	alpha = *Alpha;
	beta = *Beta;
#ಆದರೆ_ಅಂತ್ಯ

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

	ಆದರೆ (transA == 'N')
		transa = 0;
	ಆದರೆ (transA == 'T')
		transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	ಆದರೆ (transA == 'R')
		transa = 0;
	ಆದರೆ (transA == 'C')
		transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
	ಆದರೆ (transA == 'R')
		transa = 2;
	ಆದರೆ (transA == 'C')
		transa = 3;
#ಆದರೆ_ಅಂತ್ಯ

	ಆದರೆ (transB == 'N')
		transb = 0;
	ಆದರೆ (transB == 'T')
		transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	ಆದರೆ (transB == 'R')
		transb = 0;
	ಆದರೆ (transB == 'C')
		transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
	ಆದರೆ (transB == 'R')
		transb = 2;
	ಆದರೆ (transB == 'C')
		transb = 3;
#ಆದರೆ_ಅಂತ್ಯ

	ಆದರೆ (Uplo == 'U')
		uplo = 0;
	ಆದರೆ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ಆದರೆ (transa & 1) nrowa = k;
	nrowb = k;
	ಆದರೆ (transb & 1) nrowb = m;

	info = 0;

	ಆದರೆ (ldc < MAX(1, m))
		info = 13;
	ಆದರೆ (ldb < MAX(1, nrowb))
		info = 10;
	ಆದರೆ (lda < MAX(1, nrowa))
		info = 8;
	ಆದರೆ (k < 0)
		info = 5;
	ಆದರೆ (m < 0)
		info = 4;
	ಆದರೆ (transb < 0)
		info = 3;
	ಆದರೆ (transa < 0)
		info = 2;
	ಆದರೆ (uplo < 0)
		info = 1;

	ಆದರೆ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
		ಹಿಂತಿರುಗಿ;
	}
#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_UPLO Uplo,
	   ಎಣಿಕೆ CBLAS_TRANSPOSE TransA, ಎಣಿಕೆ CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint LDA,
	   ಶೂನ್ಯ *vb, blasint LDB, ಶೂನ್ಯ *vbeta, ಶೂನ್ಯ *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#ಆದರೆ_ಅಂತ್ಯ
	ಪೂರ್ಣಾಂಕ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
		ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo = 0;
		ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo = 1;

		ಆದರೆ (TransA == CblasNoTrans)
			transa = 0;
		ಆದರೆ (TransA == CblasTrans)
			transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
		ಆದರೆ (TransA == CblasConjNoTrans)
			transa = 0;
		ಆದರೆ (TransA == CblasConjTrans)
			transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
		ಆದರೆ (TransA == CblasConjNoTrans)
			transa = 2;
		ಆದರೆ (TransA == CblasConjTrans)
			transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
		ಆದರೆ (TransB == CblasNoTrans)
			transb = 0;
		ಆದರೆ (TransB == CblasTrans)
			transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
		ಆದರೆ (TransB == CblasConjNoTrans)
			transb = 0;
		ಆದರೆ (TransB == CblasConjTrans)
			transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
		ಆದರೆ (TransB == CblasConjNoTrans)
			transb = 2;
		ಆದರೆ (TransB == CblasConjTrans)
			transb = 3;
#ಆದರೆ_ಅಂತ್ಯ

		a = (ಶೂನ್ಯ *)A;
		b = (ಶೂನ್ಯ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ಆದರೆ (transa & 1) nrowa = k;
		nrowb = k;
		ಆದರೆ (transb & 1) nrowb = m;

		ಆದರೆ (ldc < MAX(1, m))
			info = 13;
		ಆದರೆ (ldb < MAX(1, nrowb))
			info = 10;
		ಆದರೆ (lda < MAX(1, nrowa))
			info = 8;
		ಆದರೆ (k < 0)
			info = 5;
		ಆದರೆ (m < 0)
			info = 4;
		ಆದರೆ (transb < 0)
			info = 3;
		ಆದರೆ (transa < 0)
			info = 2;
		ಆದರೆ (uplo < 0)
			info = 1;
	}

	ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {

		a = (ಶೂನ್ಯ *)B;
		b = (ಶೂನ್ಯ *)A;

		lda = LDB;
		ldb = LDA;

		ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo = 1;
		ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo = 0;

		ಆದರೆ (TransB == CblasNoTrans)
			transa = 0;
		ಆದರೆ (TransB == CblasTrans)
			transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
		ಆದರೆ (TransB == CblasConjNoTrans)
			transa = 0;
		ಆದರೆ (TransB == CblasConjTrans)
			transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
		ಆದರೆ (TransB == CblasConjNoTrans)
			transa = 2;
		ಆದರೆ (TransB == CblasConjTrans)
			transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
		ಆದರೆ (TransA == CblasNoTrans)
			transb = 0;
		ಆದರೆ (TransA == CblasTrans)
			transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
		ಆದರೆ (TransA == CblasConjNoTrans)
			transb = 0;
		ಆದರೆ (TransA == CblasConjTrans)
			transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
		ಆದರೆ (TransA == CblasConjNoTrans)
			transb = 2;
		ಆದರೆ (TransA == CblasConjTrans)
			transb = 3;
#ಆದರೆ_ಅಂತ್ಯ

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ಆದರೆ (transa & 1) ncola = k;
		ncolb = k;
		ಆದರೆ (transb & 1) ncolb = m;

		ಆದರೆ (ldc < MAX(1,m))
			info = 13;
		ಆದರೆ (ldb < MAX(1, ncolb))
			info = 8;
		ಆದರೆ (lda < MAX(1, ncola))
			info = 10;
		ಆದರೆ (k < 0)
			info = 5;
		ಆದರೆ (m < 0)
			info = 4;
		ಆದರೆ (transb < 0)
			info = 2;
		ಆದರೆ (transa < 0)
			info = 3;
		ಆದರೆ (uplo < 0)
			info = 1;
	}

	ಆದರೆ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
		ಹಿಂತಿರುಗಿ;
	}
#ಆದರೆ_ಅಂತ್ಯ
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿಸು GEMMT_ಕೋಷ್ಟಕ_ಸೂಚ್ಯಂಕ ((uplo << 2) | (transb << 1) | transa)
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು GEMMT_ಕೋಷ್ಟಕ_ಸೂಚ್ಯಂಕ ((uplo << 4) | (transb << 2) | transa)
#ಆದರೆ_ಅಂತ್ಯ

	blas_arg_t ವಾದಗಳು;
	FLOAT *ಬಫರ್, *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
	ದ್ವಿಗುಣ MMK;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
	ಪೂರ್ಣಾಂಕ ಬಹುಲಕ = BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
	ಪೂರ್ಣಾಂಕ ಬಹುಲಕ = BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
	ಪೂರ್ಣಾಂಕ ಬಹುಲಕ = BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
	ಪೂರ್ಣಾಂಕ ಬಹುಲಕ = BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
	ಪೂರ್ಣಾಂಕ ಬಹುಲಕ = BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
	ಪೂರ್ಣಾಂಕ ಬಹುಲಕ = BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ಇಲ್ಲದಿದ್ದರೆ
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#ಆದರೆ_ಅಂತ್ಯ
	};

	ಆದರೆ (m == 0)
		ಹಿಂತಿರುಗಿ;

	IDEBUG_ಆರಂಭ;

	FUNCTION_PROFILE_ಆರಂಭ();

	ವಾದಗಳು.m = m;
	ವಾದಗಳು.n = m;
	ವಾದಗಳು.k = k;

	ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
	ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;
	ವಾದಗಳು.c = (ಶೂನ್ಯ *)c;

	ವಾದಗಳು.lda = lda;
	ವಾದಗಳು.ldb = ldb;
	ವಾದಗಳು.ldc = ldc;

#ಆದರೆ defined(COMPLEX)
	ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)alpha;
	ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)beta;
#ಇಲ್ಲದಿದ್ದರೆ
	ವಾದಗಳು.alpha = (ಶೂನ್ಯ *)&alpha;
	ವಾದಗಳು.beta  = (ಶೂನ್ಯ *)&beta;
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
	ವಾದಗಳು.common   = ಖಾಲಿ;
	ವಾದಗಳು.nthreads = 1;
#ಆದರೆ_ಅಂತ್ಯ

	ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(0);
	ಆದರೆ (!ಬಫರ್) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
		ಹಿಂತಿರುಗಿ;
	}

	sa = (FLOAT *)((BLASLONG)ಬಫರ್ + GEMM_ಸ್ಥಾನಾಂತರ_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ಗಾತ್ರ_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ಸ್ಥಾನಾಂತರ_B);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
	MMK = (ದ್ವಿಗುಣ)(m + 1) * (ದ್ವಿಗುಣ)m * (ದ್ವಿಗುಣ)k;
	ಆದರೆ (MMK <= (SMP_ಮಿತಿ_MIN * GEMM_MULTITHREAD_ಮಿತಿ)) {
		ವಾದಗಳು.nthreads = 1;
	} ಇಲ್ಲದಿದ್ದರೆ {
		ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(3);
	}

	ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

		(gemmt[GEMMT_ಕೋಷ್ಟಕ_ಸೂಚ್ಯಂಕ])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
	} ಇಲ್ಲದಿದ್ದರೆ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		ಬಹುಲಕ |= (uplo << BLAS_UPLO_ಸರಿಸು);

		syrk_thread(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, gemmt[GEMMT_ಕೋಷ್ಟಕ_ಸೂಚ್ಯಂಕ], sa, sb, ವಾದಗಳು.nthreads);
	}
#ಆದರೆ_ಅಂತ್ಯ

	blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

	FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE, ವಾದಗಳು.m * ವಾದಗಳು.k + ವಾದಗಳು.k * ವಾದಗಳು.n + ವಾದಗಳು.m * ವಾದಗಳು.n / 2,
			     ವಾದಗಳು.m * ವಾದಗಳು.n * ವಾದಗಳು.k);

	IDEBUG_ಕೊನೆ;

	ಹಿಂತಿರುಗಿ;
}
