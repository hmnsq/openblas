/*********************************************************************/
/* Copyright 2024-2025 The OpenBLAS Project.                         */
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

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 65536.0
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SBGEMMT "

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
	IFLOAT *ਬਫਰ;
	IFLOAT *aa, *bb;
	FLOAT *cc;
	FLOAT alpha, beta;

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

	alpha = *Alpha;
	beta = *Beta;

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

	ਜੇ (transA == 'R')
		transa = 0;
	ਜੇ (transA == 'C')
		transa = 1;

	ਜੇ (transB == 'N')
		transb = 0;
	ਜੇ (transB == 'T')
		transb = 1;

	ਜੇ (transB == 'R')
		transb = 0;
	ਜੇ (transB == 'C')
		transb = 1;

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
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ਪੂਰਨ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *ਬਫਰ;

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

		ਜੇ (TransA == CblasConjNoTrans)
			transa = 0;
		ਜੇ (TransA == CblasConjTrans)
			transa = 1;

		ਜੇ (TransB == CblasNoTrans)
			transb = 0;
		ਜੇ (TransB == CblasTrans)
			transb = 1;

		ਜੇ (TransB == CblasConjNoTrans)
			transb = 0;
		ਜੇ (TransB == CblasConjTrans)
			transb = 1;

		a = (ਖਾਲੀ *)A;
		b = (ਖਾਲੀ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ਜੇ (transa & 1) nrowa = k;
		nrowb = k;
		ਜੇ (transb & 1)  nrowb = m;

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

		ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo = 0;
		ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo = 1;

		ਜੇ (TransB == CblasNoTrans)
			transa = 0;
		ਜੇ (TransB == CblasTrans)
			transa = 1;

		ਜੇ (TransB == CblasConjNoTrans)
			transa = 0;
		ਜੇ (TransB == CblasConjTrans)
			transa = 1;

		ਜੇ (TransA == CblasNoTrans)
			transb = 0;
		ਜੇ (TransA == CblasTrans)
			transb = 1;

		ਜੇ (TransA == CblasConjNoTrans)
			transb = 0;
		ਜੇ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ਜੇ (transa & 1) ncola = k;
		ncolb = k;

		ਜੇ (transb & 1) {
			ncolb = m;
		}

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
	ਪੂਰਨ ਬਫਰ_ਮਾਪ;
	blasint i, j;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
	ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ


#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
	ਸਥਿਰ ਪੂਰਨ (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ਪੂਰਨ) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#ਜੇ_ਅੰਤ
	ਪੂਰਨ (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ਜੇ (m == 0)
		ਵਾਪਸ;

	IDEBUG_ਸ਼ੁਰੂ;

	ਅਚਲ blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ਜੇ (uplo == 1) {
		ਲਈ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ਜੇ (transa & 1) {
				aa = a + lda * i;
			}
			ਜੇ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ਜੇ 0
			ਜੇ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ਨਿਰਾ, 0, ਨਿਰਾ, 0);

			ਜੇ (alpha == ZERO)
				ਜਾਰੀ;
#ਜੇ_ਅੰਤ

			IDEBUG_ਸ਼ੁਰੂ;

			ਬਫਰ_ਮਾਪ = j + k + 128 / ਆਕਾਰ(FLOAT);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ WINDOWS_ABI
			ਬਫਰ_ਮਾਪ += 160 / ਆਕਾਰ(FLOAT);
#ਜੇ_ਅੰਤ
			// for alignment
			ਬਫਰ_ਮਾਪ = (ਬਫਰ_ਮਾਪ + 3) & ~3;
			ਢੇਰ_ALLOC(ਬਫਰ_ਮਾਪ, IFLOAT, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP

			ਜੇ (1L * j * k < 2304L * GEMM_MULTITHREAD_ਹੱਦ)
				nthreads = 1;
			ਨਹੀਂਤਾਂ
				nthreads = ਅੰਕ_cpu_avail(2);

			ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

				ਜੇ (!(transa & 1))
				(gemv[(ਪੂਰਨ)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				ਨਹੀਂਤਾਂ
				(gemv[(ਪੂਰਨ)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
			} ਨਹੀਂਤਾਂ {
				ਜੇ (!(transa & 1))
				(gemv_thread[(ਪੂਰਨ)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				ਨਹੀਂਤਾਂ
				(gemv_thread[(ਪੂਰਨ)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#ਜੇ_ਅੰਤ

			ਢੇਰ_FREE(ਬਫਰ);
		}
	} ਨਹੀਂਤਾਂ {

		ਲਈ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ਜੇ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ਜੇ 0
			ਜੇ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ਨਿਰਾ, 0, ਨਿਰਾ, 0);

			ਜੇ (alpha == ZERO)
				ਜਾਰੀ;
#ਜੇ_ਅੰਤ
			IDEBUG_ਸ਼ੁਰੂ;

			ਬਫਰ_ਮਾਪ = j + k + 128 / ਆਕਾਰ(FLOAT);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ WINDOWS_ABI
			ਬਫਰ_ਮਾਪ += 160 / ਆਕਾਰ(FLOAT);
#ਜੇ_ਅੰਤ
			// for alignment
			ਬਫਰ_ਮਾਪ = (ਬਫਰ_ਮਾਪ + 3) & ~3;
			ਢੇਰ_ALLOC(ਬਫਰ_ਮਾਪ, IFLOAT, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP

			ਜੇ (1L * j * k < 2304L * GEMM_MULTITHREAD_ਹੱਦ)
				nthreads = 1;
			ਨਹੀਂਤਾਂ
				nthreads = ਅੰਕ_cpu_avail(2);

			ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

				ਜੇ (!(transa & 1))
				(gemv[(ਪੂਰਨ)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				ਨਹੀਂਤਾਂ
				(gemv[(ਪੂਰਨ)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
			} ਨਹੀਂਤਾਂ {
				ਜੇ (!(transa & 1))
				(gemv_thread[(ਪੂਰਨ)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				ਨਹੀਂਤਾਂ
				(gemv_thread[(ਪੂਰਨ)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#ਜੇ_ਅੰਤ

			ਢੇਰ_FREE(ਬਫਰ);
		}
	}

	IDEBUG_ਸਿਰਾ;

	ਵਾਪਸ;
}
