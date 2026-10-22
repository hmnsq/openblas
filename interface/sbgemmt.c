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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <stdlib.h>
#ಅನ್ತರ್ಭಾವ "common.h"

#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಮಿತಿ_MIN 65536.0
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SBGEMMT "

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
	IFLOAT *ಬಫರ್;
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

	ಆದರೆ (transA == 'N')
		transa = 0;
	ಆದರೆ (transA == 'T')
		transa = 1;

	ಆದರೆ (transA == 'R')
		transa = 0;
	ಆದರೆ (transA == 'C')
		transa = 1;

	ಆದರೆ (transB == 'N')
		transb = 0;
	ಆದರೆ (transB == 'T')
		transb = 1;

	ಆದರೆ (transB == 'R')
		transb = 0;
	ಆದರೆ (transB == 'C')
		transb = 1;

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
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ಪೂರ್ಣಾಂಕ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *ಬಫರ್;

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

		ಆದರೆ (TransA == CblasConjNoTrans)
			transa = 0;
		ಆದರೆ (TransA == CblasConjTrans)
			transa = 1;

		ಆದರೆ (TransB == CblasNoTrans)
			transb = 0;
		ಆದರೆ (TransB == CblasTrans)
			transb = 1;

		ಆದರೆ (TransB == CblasConjNoTrans)
			transb = 0;
		ಆದರೆ (TransB == CblasConjTrans)
			transb = 1;

		a = (ಶೂನ್ಯ *)A;
		b = (ಶೂನ್ಯ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ಆದರೆ (transa & 1) nrowa = k;
		nrowb = k;
		ಆದರೆ (transb & 1)  nrowb = m;

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

		ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ) uplo = 0;
		ಆದರೆ (Uplo == Cblas_ಅಧೋ) uplo = 1;

		ಆದರೆ (TransB == CblasNoTrans)
			transa = 0;
		ಆದರೆ (TransB == CblasTrans)
			transa = 1;

		ಆದರೆ (TransB == CblasConjNoTrans)
			transa = 0;
		ಆದರೆ (TransB == CblasConjTrans)
			transa = 1;

		ಆದರೆ (TransA == CblasNoTrans)
			transb = 0;
		ಆದರೆ (TransA == CblasTrans)
			transb = 1;

		ಆದರೆ (TransA == CblasConjNoTrans)
			transb = 0;
		ಆದರೆ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ಆದರೆ (transa & 1) ncola = k;
		ncolb = k;

		ಆದರೆ (transb & 1) {
			ncolb = m;
		}

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
	ಪೂರ್ಣಾಂಕ ಬಫರ್_ಗಾತ್ರ;
	blasint i, j;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
	ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ


#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
	ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ಪೂರ್ಣಾಂಕ) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#ಆದರೆ_ಅಂತ್ಯ
	ಪೂರ್ಣಾಂಕ (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ಆದರೆ (m == 0)
		ಹಿಂತಿರುಗಿ;

	IDEBUG_ಆರಂಭ;

	ಸ್ಥಿರಾಂಕ blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ಆದರೆ (uplo == 1) {
		ಗಾಗಿ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ಆದರೆ (transa & 1) {
				aa = a + lda * i;
			}
			ಆದರೆ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ಆದರೆ 0
			ಆದರೆ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ಖಾಲಿ, 0, ಖಾಲಿ, 0);

			ಆದರೆ (alpha == ZERO)
				ಮುಂದುವರಿ;
#ಆದರೆ_ಅಂತ್ಯ

			IDEBUG_ಆರಂಭ;

			ಬಫರ್_ಗಾತ್ರ = j + k + 128 / ಗಾತ್ರ(FLOAT);
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ WINDOWS_ABI
			ಬಫರ್_ಗಾತ್ರ += 160 / ಗಾತ್ರ(FLOAT);
#ಆದರೆ_ಅಂತ್ಯ
			// for alignment
			ಬಫರ್_ಗಾತ್ರ = (ಬಫರ್_ಗಾತ್ರ + 3) & ~3;
			ರಾಶಿ_ALLOC(ಬಫರ್_ಗಾತ್ರ, IFLOAT, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP

			ಆದರೆ (1L * j * k < 2304L * GEMM_MULTITHREAD_ಮಿತಿ)
				nthreads = 1;
			ಇಲ್ಲದಿದ್ದರೆ
				nthreads = ಅಂಕೆ_cpu_avail(2);

			ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

				ಆದರೆ (!(transa & 1))
				(gemv[(ಪೂರ್ಣಾಂಕ)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				ಇಲ್ಲದಿದ್ದರೆ
				(gemv[(ಪೂರ್ಣಾಂಕ)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
			} ಇಲ್ಲದಿದ್ದರೆ {
				ಆದರೆ (!(transa & 1))
				(gemv_thread[(ಪೂರ್ಣಾಂಕ)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				ಇಲ್ಲದಿದ್ದರೆ
				(gemv_thread[(ಪೂರ್ಣಾಂಕ)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#ಆದರೆ_ಅಂತ್ಯ

			ರಾಶಿ_FREE(ಬಫರ್);
		}
	} ಇಲ್ಲದಿದ್ದರೆ {

		ಗಾಗಿ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ಆದರೆ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ಆದರೆ 0
			ಆದರೆ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ಖಾಲಿ, 0, ಖಾಲಿ, 0);

			ಆದರೆ (alpha == ZERO)
				ಮುಂದುವರಿ;
#ಆದರೆ_ಅಂತ್ಯ
			IDEBUG_ಆರಂಭ;

			ಬಫರ್_ಗಾತ್ರ = j + k + 128 / ಗಾತ್ರ(FLOAT);
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ WINDOWS_ABI
			ಬಫರ್_ಗಾತ್ರ += 160 / ಗಾತ್ರ(FLOAT);
#ಆದರೆ_ಅಂತ್ಯ
			// for alignment
			ಬಫರ್_ಗಾತ್ರ = (ಬಫರ್_ಗಾತ್ರ + 3) & ~3;
			ರಾಶಿ_ALLOC(ಬಫರ್_ಗಾತ್ರ, IFLOAT, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP

			ಆದರೆ (1L * j * k < 2304L * GEMM_MULTITHREAD_ಮಿತಿ)
				nthreads = 1;
			ಇಲ್ಲದಿದ್ದರೆ
				nthreads = ಅಂಕೆ_cpu_avail(2);

			ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

				ಆದರೆ (!(transa & 1))
				(gemv[(ಪೂರ್ಣಾಂಕ)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				ಇಲ್ಲದಿದ್ದರೆ
				(gemv[(ಪೂರ್ಣಾಂಕ)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
			} ಇಲ್ಲದಿದ್ದರೆ {
				ಆದರೆ (!(transa & 1))
				(gemv_thread[(ಪೂರ್ಣಾಂಕ)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				ಇಲ್ಲದಿದ್ದರೆ
				(gemv_thread[(ಪೂರ್ಣಾಂಕ)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#ಆದರೆ_ಅಂತ್ಯ

			ರಾಶಿ_FREE(ಬಫರ್);
		}
	}

	IDEBUG_ಕೊನೆ;

	ಹಿಂತಿರುಗಿ;
}
