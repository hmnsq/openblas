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

void NAME(var_na *UPLO, var_na *TRANSA, var_na *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ចំនួនគត់ transa, transb, uplo;
	blasint info;

	var_na transA, transB, Uplo;
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

	ប្រសិនបើ (transA == 'N')
		transa = 0;
	ប្រសិនបើ (transA == 'T')
		transa = 1;

	ប្រសិនបើ (transA == 'R')
		transa = 0;
	ប្រសិនបើ (transA == 'C')
		transa = 1;

	ប្រសិនបើ (transB == 'N')
		transb = 0;
	ប្រសិនបើ (transB == 'T')
		transb = 1;

	ប្រសិនបើ (transB == 'R')
		transb = 0;
	ប្រសិនបើ (transB == 'C')
		transb = 1;

	ប្រសិនបើ (Uplo == 'U')
		uplo = 0;
	ប្រសិនបើ (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	ប្រសិនបើ (transa & 1) nrowa = k;
	nrowb = k;
	ប្រសិនបើ (transb & 1) nrowb = m;

	info = 0;

	ប្រសិនបើ (ldc < MAX(1, m))
		info = 13;
	ប្រសិនបើ (ldb < MAX(1, nrowb))
		info = 10;
	ប្រសិនបើ (lda < MAX(1, nrowa))
		info = 8;
	ប្រសិនបើ (k < 0)
		info = 5;
	ប្រសិនបើ (m < 0)
		info = 4;
	ប្រសិនបើ (transb < 0)
		info = 3;
	ប្រសិនបើ (transa < 0)
		info = 2;
	ប្រសិនបើ (uplo < 0)
		info = 1;

	ប្រសិនបើ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ត្រឡប់;
	}
#បើមិនដូច្នេះ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ចំនួនគត់ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ប្រសិនបើ (order == CblasColMajor) {
		ប្រសិនបើ (Uplo == CblasUpper) uplo = 0;
		ប្រសិនបើ (Uplo == CblasLower) uplo = 1;

		ប្រសិនបើ (TransA == CblasNoTrans)
			transa = 0;
		ប្រសិនបើ (TransA == CblasTrans)
			transa = 1;

		ប្រសិនបើ (TransA == CblasConjNoTrans)
			transa = 0;
		ប្រសិនបើ (TransA == CblasConjTrans)
			transa = 1;

		ប្រសិនបើ (TransB == CblasNoTrans)
			transb = 0;
		ប្រសិនបើ (TransB == CblasTrans)
			transb = 1;

		ប្រសិនបើ (TransB == CblasConjNoTrans)
			transb = 0;
		ប្រសិនបើ (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ប្រសិនបើ (transa & 1) nrowa = k;
		nrowb = k;
		ប្រសិនបើ (transb & 1)  nrowb = m;

		ប្រសិនបើ (ldc < MAX(1, m))
			info = 13;
		ប្រសិនបើ (ldb < MAX(1, nrowb))
			info = 10;
		ប្រសិនបើ (lda < MAX(1, nrowa))
			info = 8;
		ប្រសិនបើ (k < 0)
			info = 5;
		ប្រសិនបើ (m < 0)
			info = 4;
		ប្រសិនបើ (transb < 0)
			info = 3;
		ប្រសិនបើ (transa < 0)
			info = 2;
		ប្រសិនបើ (uplo < 0)
			info = 1;
	}

	ប្រសិនបើ (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		ប្រសិនបើ (Uplo == CblasUpper) uplo = 0;
		ប្រសិនបើ (Uplo == CblasLower) uplo = 1;

		ប្រសិនបើ (TransB == CblasNoTrans)
			transa = 0;
		ប្រសិនបើ (TransB == CblasTrans)
			transa = 1;

		ប្រសិនបើ (TransB == CblasConjNoTrans)
			transa = 0;
		ប្រសិនបើ (TransB == CblasConjTrans)
			transa = 1;

		ប្រសិនបើ (TransA == CblasNoTrans)
			transb = 0;
		ប្រសិនបើ (TransA == CblasTrans)
			transb = 1;

		ប្រសិនបើ (TransA == CblasConjNoTrans)
			transb = 0;
		ប្រសិនបើ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ប្រសិនបើ (transa & 1) ncola = k;
		ncolb = k;

		ប្រសិនបើ (transb & 1) {
			ncolb = m;
		}

		ប្រសិនបើ (ldc < MAX(1,m))
			info = 13;
		ប្រសិនបើ (ldb < MAX(1, ncolb))
			info = 8;
		ប្រសិនបើ (lda < MAX(1, ncola))
			info = 10;
		ប្រសិនបើ (k < 0)
			info = 5;
		ប្រសិនបើ (m < 0)
			info = 4;
		ប្រសិនបើ (transb < 0)
			info = 2;
		ប្រសិនបើ (transa < 0)
			info = 3;
		ប្រសិនបើ (uplo < 0)
			info = 1;
	}

	ប្រសិនបើ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ត្រឡប់;
	}

#endif
	ចំនួនគត់ buffer_size;
	blasint i, j;

#ifdef SMP
	ចំនួនគត់ nthreads;
#endif


#ifdef SMP
	static ចំនួនគត់ (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ចំនួនគត់) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	ចំនួនគត់ (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ប្រសិនបើ (m == 0)
		ត្រឡប់;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ប្រសិនបើ (uplo == 1) {
		សម្រាប់ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ប្រសិនបើ (transa & 1) {
				aa = a + lda * i;
			}
			ប្រសិនបើ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ប្រសិនបើ 0
			ប្រសិនបើ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			ប្រសិនបើ (alpha == ZERO)
				continue;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / parimaa_na(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / parimaa_na(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			STACK_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ប្រសិនបើ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			បើមិនដូច្នេះ
				nthreads = num_cpu_avail(2);

			ប្រសិនបើ (nthreads == 1) {
#endif

				ប្រសិនបើ (!(transa & 1))
				(gemv[(ចំនួនគត់)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				បើមិនដូច្នេះ
				(gemv[(ចំនួនគត់)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} បើមិនដូច្នេះ {
				ប្រសិនបើ (!(transa & 1))
				(gemv_thread[(ចំនួនគត់)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				បើមិនដូច្នេះ
				(gemv_thread[(ចំនួនគត់)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} បើមិនដូច្នេះ {

		សម្រាប់ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ប្រសិនបើ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ប្រសិនបើ 0
			ប្រសិនបើ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			ប្រសិនបើ (alpha == ZERO)
				continue;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / parimaa_na(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / parimaa_na(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			STACK_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ប្រសិនបើ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			បើមិនដូច្នេះ
				nthreads = num_cpu_avail(2);

			ប្រសិនបើ (nthreads == 1) {
#endif

				ប្រសិនបើ (!(transa & 1))
				(gemv[(ចំនួនគត់)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				បើមិនដូច្នេះ
				(gemv[(ចំនួនគត់)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} បើមិនដូច្នេះ {
				ប្រសិនបើ (!(transa & 1))
				(gemv_thread[(ចំនួនគត់)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				បើមិនដូច្នេះ
				(gemv_thread[(ចំនួនគត់)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	ត្រឡប់;
}
