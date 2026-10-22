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

#include <stdio.h>
#include <stdlib.h>
#include "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, char *TRANSA, char *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	цяло transa, transb, uplo;
	blasint info;

	char transA, transB, Uplo;
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

	ако (transA == 'N')
		transa = 0;
	ако (transA == 'T')
		transa = 1;

	ако (transA == 'R')
		transa = 0;
	ако (transA == 'C')
		transa = 1;

	ако (transB == 'N')
		transb = 0;
	ако (transB == 'T')
		transb = 1;

	ако (transB == 'R')
		transb = 0;
	ако (transB == 'C')
		transb = 1;

	ако (Uplo == 'U')
		uplo = 0;
	ако (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	ако (transa & 1) nrowa = k;
	nrowb = k;
	ако (transb & 1) nrowb = m;

	info = 0;

	ако (ldc < MAX(1, m))
		info = 13;
	ако (ldb < MAX(1, nrowb))
		info = 10;
	ако (lda < MAX(1, nrowa))
		info = 8;
	ако (k < 0)
		info = 5;
	ако (m < 0)
		info = 4;
	ако (transb < 0)
		info = 3;
	ако (transa < 0)
		info = 2;
	ако (uplo < 0)
		info = 1;

	ако (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, sizeof(ERROR_NAME));
		върни;
	}
#иначе

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo,
	   enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	цяло transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ако (order == CblasColMajor) {
		ако (Uplo == CblasUpper) uplo = 0;
		ако (Uplo == CblasLower) uplo = 1;

		ако (TransA == CblasNoTrans)
			transa = 0;
		ако (TransA == CblasTrans)
			transa = 1;

		ако (TransA == CblasConjNoTrans)
			transa = 0;
		ако (TransA == CblasConjTrans)
			transa = 1;

		ако (TransB == CblasNoTrans)
			transb = 0;
		ако (TransB == CblasTrans)
			transb = 1;

		ако (TransB == CblasConjNoTrans)
			transb = 0;
		ако (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ако (transa & 1) nrowa = k;
		nrowb = k;
		ако (transb & 1)  nrowb = m;

		ако (ldc < MAX(1, m))
			info = 13;
		ако (ldb < MAX(1, nrowb))
			info = 10;
		ако (lda < MAX(1, nrowa))
			info = 8;
		ако (k < 0)
			info = 5;
		ако (m < 0)
			info = 4;
		ако (transb < 0)
			info = 3;
		ако (transa < 0)
			info = 2;
		ако (uplo < 0)
			info = 1;
	}

	ако (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		ако (Uplo == CblasUpper) uplo = 0;
		ако (Uplo == CblasLower) uplo = 1;

		ако (TransB == CblasNoTrans)
			transa = 0;
		ако (TransB == CblasTrans)
			transa = 1;

		ако (TransB == CblasConjNoTrans)
			transa = 0;
		ако (TransB == CblasConjTrans)
			transa = 1;

		ако (TransA == CblasNoTrans)
			transb = 0;
		ако (TransA == CblasTrans)
			transb = 1;

		ако (TransA == CblasConjNoTrans)
			transb = 0;
		ако (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ако (transa & 1) ncola = k;
		ncolb = k;

		ако (transb & 1) {
			ncolb = m;
		}

		ако (ldc < MAX(1,m))
			info = 13;
		ако (ldb < MAX(1, ncolb))
			info = 8;
		ако (lda < MAX(1, ncola))
			info = 10;
		ако (k < 0)
			info = 5;
		ако (m < 0)
			info = 4;
		ако (transb < 0)
			info = 2;
		ако (transa < 0)
			info = 3;
		ако (uplo < 0)
			info = 1;
	}

	ако (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, sizeof(ERROR_NAME));
		върни;
	}

#endif
	цяло buffer_size;
	blasint i, j;

#ifdef SMP
	цяло nthreads;
#endif


#ifdef SMP
	static цяло (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, цяло) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	цяло (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ако (m == 0)
		върни;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ако (uplo == 1) {
		за (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ако (transa & 1) {
				aa = a + lda * i;
			}
			ако (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ако 0
			ако (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			ако (alpha == ZERO)
				continue;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / sizeof(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / sizeof(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			STACK_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ако (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			иначе
				nthreads = num_cpu_avail(2);

			ако (nthreads == 1) {
#endif

				ако (!(transa & 1))
				(gemv[(цяло)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				иначе
				(gemv[(цяло)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} иначе {
				ако (!(transa & 1))
				(gemv_thread[(цяло)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				иначе
				(gemv_thread[(цяло)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} иначе {

		за (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ако (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ако 0
			ако (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			ако (alpha == ZERO)
				continue;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / sizeof(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / sizeof(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			STACK_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ако (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			иначе
				nthreads = num_cpu_avail(2);

			ако (nthreads == 1) {
#endif

				ако (!(transa & 1))
				(gemv[(цяло)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				иначе
				(gemv[(цяло)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} иначе {
				ако (!(transa & 1))
				(gemv_thread[(цяло)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				иначе
				(gemv_thread[(цяло)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	върни;
}
