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
	ціле transa, transb, uplo;
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

	якщо (transA == 'N')
		transa = 0;
	якщо (transA == 'T')
		transa = 1;

	якщо (transA == 'R')
		transa = 0;
	якщо (transA == 'C')
		transa = 1;

	якщо (transB == 'N')
		transb = 0;
	якщо (transB == 'T')
		transb = 1;

	якщо (transB == 'R')
		transb = 0;
	якщо (transB == 'C')
		transb = 1;

	якщо (Uplo == 'U')
		uplo = 0;
	якщо (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	якщо (transa & 1) nrowa = k;
	nrowb = k;
	якщо (transb & 1) nrowb = m;

	info = 0;

	якщо (ldc < MAX(1, m))
		info = 13;
	якщо (ldb < MAX(1, nrowb))
		info = 10;
	якщо (lda < MAX(1, nrowa))
		info = 8;
	якщо (k < 0)
		info = 5;
	якщо (m < 0)
		info = 4;
	якщо (transb < 0)
		info = 3;
	якщо (transa < 0)
		info = 2;
	якщо (uplo < 0)
		info = 1;

	якщо (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, sizeof(ERROR_NAME));
		повернути;
	}
#інакше

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo,
	   enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ціле transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	якщо (order == CblasColMajor) {
		якщо (Uplo == CblasUpper) uplo = 0;
		якщо (Uplo == CblasLower) uplo = 1;

		якщо (TransA == CblasNoTrans)
			transa = 0;
		якщо (TransA == CblasTrans)
			transa = 1;

		якщо (TransA == CblasConjNoTrans)
			transa = 0;
		якщо (TransA == CblasConjTrans)
			transa = 1;

		якщо (TransB == CblasNoTrans)
			transb = 0;
		якщо (TransB == CblasTrans)
			transb = 1;

		якщо (TransB == CblasConjNoTrans)
			transb = 0;
		якщо (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		якщо (transa & 1) nrowa = k;
		nrowb = k;
		якщо (transb & 1)  nrowb = m;

		якщо (ldc < MAX(1, m))
			info = 13;
		якщо (ldb < MAX(1, nrowb))
			info = 10;
		якщо (lda < MAX(1, nrowa))
			info = 8;
		якщо (k < 0)
			info = 5;
		якщо (m < 0)
			info = 4;
		якщо (transb < 0)
			info = 3;
		якщо (transa < 0)
			info = 2;
		якщо (uplo < 0)
			info = 1;
	}

	якщо (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		якщо (Uplo == CblasUpper) uplo = 0;
		якщо (Uplo == CblasLower) uplo = 1;

		якщо (TransB == CblasNoTrans)
			transa = 0;
		якщо (TransB == CblasTrans)
			transa = 1;

		якщо (TransB == CblasConjNoTrans)
			transa = 0;
		якщо (TransB == CblasConjTrans)
			transa = 1;

		якщо (TransA == CblasNoTrans)
			transb = 0;
		якщо (TransA == CblasTrans)
			transb = 1;

		якщо (TransA == CblasConjNoTrans)
			transb = 0;
		якщо (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		якщо (transa & 1) ncola = k;
		ncolb = k;

		якщо (transb & 1) {
			ncolb = m;
		}

		якщо (ldc < MAX(1,m))
			info = 13;
		якщо (ldb < MAX(1, ncolb))
			info = 8;
		якщо (lda < MAX(1, ncola))
			info = 10;
		якщо (k < 0)
			info = 5;
		якщо (m < 0)
			info = 4;
		якщо (transb < 0)
			info = 2;
		якщо (transa < 0)
			info = 3;
		якщо (uplo < 0)
			info = 1;
	}

	якщо (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, sizeof(ERROR_NAME));
		повернути;
	}

#endif
	ціле buffer_size;
	blasint i, j;

#ifdef SMP
	ціле nthreads;
#endif


#ifdef SMP
	static ціле (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ціле) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	ціле (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	якщо (m == 0)
		повернути;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	якщо (uplo == 1) {
		для (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			якщо (transa & 1) {
				aa = a + lda * i;
			}
			якщо (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#якщо 0
			якщо (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			якщо (alpha == ZERO)
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

			якщо (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			інакше
				nthreads = num_cpu_avail(2);

			якщо (nthreads == 1) {
#endif

				якщо (!(transa & 1))
				(gemv[(ціле)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				інакше
				(gemv[(ціле)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} інакше {
				якщо (!(transa & 1))
				(gemv_thread[(ціле)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				інакше
				(gemv_thread[(ціле)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} інакше {

		для (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			якщо (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#якщо 0
			якщо (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			якщо (alpha == ZERO)
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

			якщо (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			інакше
				nthreads = num_cpu_avail(2);

			якщо (nthreads == 1) {
#endif

				якщо (!(transa & 1))
				(gemv[(ціле)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				інакше
				(gemv[(ціле)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} інакше {
				якщо (!(transa & 1))
				(gemv_thread[(ціле)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				інакше
				(gemv_thread[(ціле)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	повернути;
}
