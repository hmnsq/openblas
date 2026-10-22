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
	целое transa, transb, uplo;
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

	если (transA == 'N')
		transa = 0;
	если (transA == 'T')
		transa = 1;

	если (transA == 'R')
		transa = 0;
	если (transA == 'C')
		transa = 1;

	если (transB == 'N')
		transb = 0;
	если (transB == 'T')
		transb = 1;

	если (transB == 'R')
		transb = 0;
	если (transB == 'C')
		transb = 1;

	если (Uplo == 'U')
		uplo = 0;
	если (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	если (transa & 1) nrowa = k;
	nrowb = k;
	если (transb & 1) nrowb = m;

	info = 0;

	если (ldc < MAX(1, m))
		info = 13;
	если (ldb < MAX(1, nrowb))
		info = 10;
	если (lda < MAX(1, nrowa))
		info = 8;
	если (k < 0)
		info = 5;
	если (m < 0)
		info = 4;
	если (transb < 0)
		info = 3;
	если (transa < 0)
		info = 2;
	если (uplo < 0)
		info = 1;

	если (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, sizeof(ERROR_NAME));
		вернуть;
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

	целое transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	если (order == CblasColMajor) {
		если (Uplo == CblasUpper) uplo = 0;
		если (Uplo == CblasLower) uplo = 1;

		если (TransA == CblasNoTrans)
			transa = 0;
		если (TransA == CblasTrans)
			transa = 1;

		если (TransA == CblasConjNoTrans)
			transa = 0;
		если (TransA == CblasConjTrans)
			transa = 1;

		если (TransB == CblasNoTrans)
			transb = 0;
		если (TransB == CblasTrans)
			transb = 1;

		если (TransB == CblasConjNoTrans)
			transb = 0;
		если (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		если (transa & 1) nrowa = k;
		nrowb = k;
		если (transb & 1)  nrowb = m;

		если (ldc < MAX(1, m))
			info = 13;
		если (ldb < MAX(1, nrowb))
			info = 10;
		если (lda < MAX(1, nrowa))
			info = 8;
		если (k < 0)
			info = 5;
		если (m < 0)
			info = 4;
		если (transb < 0)
			info = 3;
		если (transa < 0)
			info = 2;
		если (uplo < 0)
			info = 1;
	}

	если (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		если (Uplo == CblasUpper) uplo = 0;
		если (Uplo == CblasLower) uplo = 1;

		если (TransB == CblasNoTrans)
			transa = 0;
		если (TransB == CblasTrans)
			transa = 1;

		если (TransB == CblasConjNoTrans)
			transa = 0;
		если (TransB == CblasConjTrans)
			transa = 1;

		если (TransA == CblasNoTrans)
			transb = 0;
		если (TransA == CblasTrans)
			transb = 1;

		если (TransA == CblasConjNoTrans)
			transb = 0;
		если (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		если (transa & 1) ncola = k;
		ncolb = k;

		если (transb & 1) {
			ncolb = m;
		}

		если (ldc < MAX(1,m))
			info = 13;
		если (ldb < MAX(1, ncolb))
			info = 8;
		если (lda < MAX(1, ncola))
			info = 10;
		если (k < 0)
			info = 5;
		если (m < 0)
			info = 4;
		если (transb < 0)
			info = 2;
		если (transa < 0)
			info = 3;
		если (uplo < 0)
			info = 1;
	}

	если (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, sizeof(ERROR_NAME));
		вернуть;
	}

#endif
	целое buffer_size;
	blasint i, j;

#ifdef SMP
	целое nthreads;
#endif


#ifdef SMP
	static целое (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, целое) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	целое (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	если (m == 0)
		вернуть;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	если (uplo == 1) {
		для (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			если (transa & 1) {
				aa = a + lda * i;
			}
			если (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#если 0
			если (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			если (alpha == ZERO)
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

			если (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			иначе
				nthreads = num_cpu_avail(2);

			если (nthreads == 1) {
#endif

				если (!(transa & 1))
				(gemv[(целое)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				иначе
				(gemv[(целое)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} иначе {
				если (!(transa & 1))
				(gemv_thread[(целое)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				иначе
				(gemv_thread[(целое)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} иначе {

		для (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			если (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#если 0
			если (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			если (alpha == ZERO)
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

			если (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			иначе
				nthreads = num_cpu_avail(2);

			если (nthreads == 1) {
#endif

				если (!(transa & 1))
				(gemv[(целое)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				иначе
				(gemv[(целое)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} иначе {
				если (!(transa & 1))
				(gemv_thread[(целое)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				иначе
				(gemv_thread[(целое)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	вернуть;
}
