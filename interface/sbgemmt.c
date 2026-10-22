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
	ཧྲིལ transa, transb, uplo;
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

	གལ་ཏེ (transA == 'N')
		transa = 0;
	གལ་ཏེ (transA == 'T')
		transa = 1;

	གལ་ཏེ (transA == 'R')
		transa = 0;
	གལ་ཏེ (transA == 'C')
		transa = 1;

	གལ་ཏེ (transB == 'N')
		transb = 0;
	གལ་ཏེ (transB == 'T')
		transb = 1;

	གལ་ཏེ (transB == 'R')
		transb = 0;
	གལ་ཏེ (transB == 'C')
		transb = 1;

	གལ་ཏེ (Uplo == 'U')
		uplo = 0;
	གལ་ཏེ (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	གལ་ཏེ (transa & 1) nrowa = k;
	nrowb = k;
	གལ་ཏེ (transb & 1) nrowb = m;

	info = 0;

	གལ་ཏེ (ldc < MAX(1, m))
		info = 13;
	གལ་ཏེ (ldb < MAX(1, nrowb))
		info = 10;
	གལ་ཏེ (lda < MAX(1, nrowa))
		info = 8;
	གལ་ཏེ (k < 0)
		info = 5;
	གལ་ཏེ (m < 0)
		info = 4;
	གལ་ཏེ (transb < 0)
		info = 3;
	གལ་ཏེ (transa < 0)
		info = 2;
	གལ་ཏེ (uplo < 0)
		info = 1;

	གལ་ཏེ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ཕྱིར;
	}
#དེ་མིན

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ཧྲིལ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	གལ་ཏེ (order == CblasColMajor) {
		གལ་ཏེ (Uplo == CblasUpper) uplo = 0;
		གལ་ཏེ (Uplo == CblasLower) uplo = 1;

		གལ་ཏེ (TransA == CblasNoTrans)
			transa = 0;
		གལ་ཏེ (TransA == CblasTrans)
			transa = 1;

		གལ་ཏེ (TransA == CblasConjNoTrans)
			transa = 0;
		གལ་ཏེ (TransA == CblasConjTrans)
			transa = 1;

		གལ་ཏེ (TransB == CblasNoTrans)
			transb = 0;
		གལ་ཏེ (TransB == CblasTrans)
			transb = 1;

		གལ་ཏེ (TransB == CblasConjNoTrans)
			transb = 0;
		གལ་ཏེ (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		གལ་ཏེ (transa & 1) nrowa = k;
		nrowb = k;
		གལ་ཏེ (transb & 1)  nrowb = m;

		གལ་ཏེ (ldc < MAX(1, m))
			info = 13;
		གལ་ཏེ (ldb < MAX(1, nrowb))
			info = 10;
		གལ་ཏེ (lda < MAX(1, nrowa))
			info = 8;
		གལ་ཏེ (k < 0)
			info = 5;
		གལ་ཏེ (m < 0)
			info = 4;
		གལ་ཏེ (transb < 0)
			info = 3;
		གལ་ཏེ (transa < 0)
			info = 2;
		གལ་ཏེ (uplo < 0)
			info = 1;
	}

	གལ་ཏེ (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		གལ་ཏེ (Uplo == CblasUpper) uplo = 0;
		གལ་ཏེ (Uplo == CblasLower) uplo = 1;

		གལ་ཏེ (TransB == CblasNoTrans)
			transa = 0;
		གལ་ཏེ (TransB == CblasTrans)
			transa = 1;

		གལ་ཏེ (TransB == CblasConjNoTrans)
			transa = 0;
		གལ་ཏེ (TransB == CblasConjTrans)
			transa = 1;

		གལ་ཏེ (TransA == CblasNoTrans)
			transb = 0;
		གལ་ཏེ (TransA == CblasTrans)
			transb = 1;

		གལ་ཏེ (TransA == CblasConjNoTrans)
			transb = 0;
		གལ་ཏེ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		གལ་ཏེ (transa & 1) ncola = k;
		ncolb = k;

		གལ་ཏེ (transb & 1) {
			ncolb = m;
		}

		གལ་ཏེ (ldc < MAX(1,m))
			info = 13;
		གལ་ཏེ (ldb < MAX(1, ncolb))
			info = 8;
		གལ་ཏེ (lda < MAX(1, ncola))
			info = 10;
		གལ་ཏེ (k < 0)
			info = 5;
		གལ་ཏེ (m < 0)
			info = 4;
		གལ་ཏེ (transb < 0)
			info = 2;
		གལ་ཏེ (transa < 0)
			info = 3;
		གལ་ཏེ (uplo < 0)
			info = 1;
	}

	གལ་ཏེ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ཕྱིར;
	}

#endif
	ཧྲིལ buffer_size;
	blasint i, j;

#ifdef SMP
	ཧྲིལ nthreads;
#endif


#ifdef SMP
	static ཧྲིལ (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ཧྲིལ) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	ཧྲིལ (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	གལ་ཏེ (m == 0)
		ཕྱིར;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	གལ་ཏེ (uplo == 1) {
		ཆེད་དུ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			གལ་ཏེ (transa & 1) {
				aa = a + lda * i;
			}
			གལ་ཏེ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#གལ་ཏེ 0
			གལ་ཏེ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			གལ་ཏེ (alpha == ZERO)
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

			གལ་ཏེ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			དེ་མིན
				nthreads = num_cpu_avail(2);

			གལ་ཏེ (nthreads == 1) {
#endif

				གལ་ཏེ (!(transa & 1))
				(gemv[(ཧྲིལ)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				དེ་མིན
				(gemv[(ཧྲིལ)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} དེ་མིན {
				གལ་ཏེ (!(transa & 1))
				(gemv_thread[(ཧྲིལ)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				དེ་མིན
				(gemv_thread[(ཧྲིལ)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} དེ་མིན {

		ཆེད་དུ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			གལ་ཏེ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#གལ་ཏེ 0
			གལ་ཏེ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			གལ་ཏེ (alpha == ZERO)
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

			གལ་ཏེ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			དེ་མིན
				nthreads = num_cpu_avail(2);

			གལ་ཏེ (nthreads == 1) {
#endif

				གལ་ཏེ (!(transa & 1))
				(gemv[(ཧྲིལ)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				དེ་མིན
				(gemv[(ཧྲིལ)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} དེ་མིན {
				གལ་ཏེ (!(transa & 1))
				(gemv_thread[(ཧྲིལ)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				དེ་མིན
				(gemv_thread[(ཧྲིལ)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	ཕྱིར;
}
