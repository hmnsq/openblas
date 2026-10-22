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

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

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
	IFLOAT *buffer;
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

	ಯದಿ (transA == 'N')
		transa = 0;
	ಯದಿ (transA == 'T')
		transa = 1;

	ಯದಿ (transA == 'R')
		transa = 0;
	ಯದಿ (transA == 'C')
		transa = 1;

	ಯದಿ (transB == 'N')
		transb = 0;
	ಯದಿ (transB == 'T')
		transb = 1;

	ಯದಿ (transB == 'R')
		transb = 0;
	ಯದಿ (transB == 'C')
		transb = 1;

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
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ಪೂರ್ಣಾಙ್ಕ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

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

		ಯದಿ (TransA == CblasConjNoTrans)
			transa = 0;
		ಯದಿ (TransA == CblasConjTrans)
			transa = 1;

		ಯದಿ (TransB == CblasNoTrans)
			transb = 0;
		ಯದಿ (TransB == CblasTrans)
			transb = 1;

		ಯದಿ (TransB == CblasConjNoTrans)
			transb = 0;
		ಯದಿ (TransB == CblasConjTrans)
			transb = 1;

		a = (ಶೂನ್ಯ *)A;
		b = (ಶೂನ್ಯ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ಯದಿ (transa & 1) nrowa = k;
		nrowb = k;
		ಯದಿ (transb & 1)  nrowb = m;

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

		ಯದಿ (Uplo == CblasUpper) uplo = 0;
		ಯದಿ (Uplo == CblasLower) uplo = 1;

		ಯದಿ (TransB == CblasNoTrans)
			transa = 0;
		ಯದಿ (TransB == CblasTrans)
			transa = 1;

		ಯದಿ (TransB == CblasConjNoTrans)
			transa = 0;
		ಯದಿ (TransB == CblasConjTrans)
			transa = 1;

		ಯದಿ (TransA == CblasNoTrans)
			transb = 0;
		ಯದಿ (TransA == CblasTrans)
			transb = 1;

		ಯದಿ (TransA == CblasConjNoTrans)
			transb = 0;
		ಯದಿ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ಯದಿ (transa & 1) ncola = k;
		ncolb = k;

		ಯದಿ (transb & 1) {
			ncolb = m;
		}

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
	ಪೂರ್ಣಾಙ್ಕ buffer_size;
	blasint i, j;

#ifdef SMP
	ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif


#ifdef SMP
	ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ಪೂರ್ಣಾಙ್ಕ) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	ಪೂರ್ಣಾಙ್ಕ (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ಯದಿ (m == 0)
		ಪ್ರತ್ಯಾಗಮ;

	IDEBUG_START;

	ಧ್ರುವ blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ಯದಿ (uplo == 1) {
		ಪ್ರತಿ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ಯದಿ (transa & 1) {
				aa = a + lda * i;
			}
			ಯದಿ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ಯದಿ 0
			ಯದಿ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

			ಯದಿ (alpha == ZERO)
				ಅನುವರ್ತ;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / ಪರಿಮಾಣ(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / ಪರಿಮಾಣ(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			ರಾಶಿ_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ಯದಿ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			ಅನ್ಯಥಾ
				nthreads = ಅಂಕೆ_cpu_avail(2);

			ಯದಿ (nthreads == 1) {
#endif

				ಯದಿ (!(transa & 1))
				(gemv[(ಪೂರ್ಣಾಙ್ಕ)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				ಅನ್ಯಥಾ
				(gemv[(ಪೂರ್ಣಾಙ್ಕ)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} ಅನ್ಯಥಾ {
				ಯದಿ (!(transa & 1))
				(gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				ಅನ್ಯಥಾ
				(gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			ರಾಶಿ_FREE(buffer);
		}
	} ಅನ್ಯಥಾ {

		ಪ್ರತಿ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ಯದಿ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ಯದಿ 0
			ಯದಿ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

			ಯದಿ (alpha == ZERO)
				ಅನುವರ್ತ;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / ಪರಿಮಾಣ(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / ಪರಿಮಾಣ(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			ರಾಶಿ_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ಯದಿ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			ಅನ್ಯಥಾ
				nthreads = ಅಂಕೆ_cpu_avail(2);

			ಯದಿ (nthreads == 1) {
#endif

				ಯದಿ (!(transa & 1))
				(gemv[(ಪೂರ್ಣಾಙ್ಕ)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				ಅನ್ಯಥಾ
				(gemv[(ಪೂರ್ಣಾಙ್ಕ)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} ಅನ್ಯಥಾ {
				ಯದಿ (!(transa & 1))
				(gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				ಅನ್ಯಥಾ
				(gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			ರಾಶಿ_FREE(buffer);
		}
	}

	IDEBUG_END;

	ಪ್ರತ್ಯಾಗಮ;
}
