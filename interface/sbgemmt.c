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
	จำนวนเต็ม transa, transb, uplo;
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

	ถ้า (transA == 'N')
		transa = 0;
	ถ้า (transA == 'T')
		transa = 1;

	ถ้า (transA == 'R')
		transa = 0;
	ถ้า (transA == 'C')
		transa = 1;

	ถ้า (transB == 'N')
		transb = 0;
	ถ้า (transB == 'T')
		transb = 1;

	ถ้า (transB == 'R')
		transb = 0;
	ถ้า (transB == 'C')
		transb = 1;

	ถ้า (Uplo == 'U')
		uplo = 0;
	ถ้า (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	ถ้า (transa & 1) nrowa = k;
	nrowb = k;
	ถ้า (transb & 1) nrowb = m;

	info = 0;

	ถ้า (ldc < MAX(1, m))
		info = 13;
	ถ้า (ldb < MAX(1, nrowb))
		info = 10;
	ถ้า (lda < MAX(1, nrowa))
		info = 8;
	ถ้า (k < 0)
		info = 5;
	ถ้า (m < 0)
		info = 4;
	ถ้า (transb < 0)
		info = 3;
	ถ้า (transa < 0)
		info = 2;
	ถ้า (uplo < 0)
		info = 1;

	ถ้า (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		คืน;
	}
#มิฉะนั้น

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	จำนวนเต็ม transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ถ้า (order == CblasColMajor) {
		ถ้า (Uplo == CblasUpper) uplo = 0;
		ถ้า (Uplo == CblasLower) uplo = 1;

		ถ้า (TransA == CblasNoTrans)
			transa = 0;
		ถ้า (TransA == CblasTrans)
			transa = 1;

		ถ้า (TransA == CblasConjNoTrans)
			transa = 0;
		ถ้า (TransA == CblasConjTrans)
			transa = 1;

		ถ้า (TransB == CblasNoTrans)
			transb = 0;
		ถ้า (TransB == CblasTrans)
			transb = 1;

		ถ้า (TransB == CblasConjNoTrans)
			transb = 0;
		ถ้า (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ถ้า (transa & 1) nrowa = k;
		nrowb = k;
		ถ้า (transb & 1)  nrowb = m;

		ถ้า (ldc < MAX(1, m))
			info = 13;
		ถ้า (ldb < MAX(1, nrowb))
			info = 10;
		ถ้า (lda < MAX(1, nrowa))
			info = 8;
		ถ้า (k < 0)
			info = 5;
		ถ้า (m < 0)
			info = 4;
		ถ้า (transb < 0)
			info = 3;
		ถ้า (transa < 0)
			info = 2;
		ถ้า (uplo < 0)
			info = 1;
	}

	ถ้า (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		ถ้า (Uplo == CblasUpper) uplo = 0;
		ถ้า (Uplo == CblasLower) uplo = 1;

		ถ้า (TransB == CblasNoTrans)
			transa = 0;
		ถ้า (TransB == CblasTrans)
			transa = 1;

		ถ้า (TransB == CblasConjNoTrans)
			transa = 0;
		ถ้า (TransB == CblasConjTrans)
			transa = 1;

		ถ้า (TransA == CblasNoTrans)
			transb = 0;
		ถ้า (TransA == CblasTrans)
			transb = 1;

		ถ้า (TransA == CblasConjNoTrans)
			transb = 0;
		ถ้า (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ถ้า (transa & 1) ncola = k;
		ncolb = k;

		ถ้า (transb & 1) {
			ncolb = m;
		}

		ถ้า (ldc < MAX(1,m))
			info = 13;
		ถ้า (ldb < MAX(1, ncolb))
			info = 8;
		ถ้า (lda < MAX(1, ncola))
			info = 10;
		ถ้า (k < 0)
			info = 5;
		ถ้า (m < 0)
			info = 4;
		ถ้า (transb < 0)
			info = 2;
		ถ้า (transa < 0)
			info = 3;
		ถ้า (uplo < 0)
			info = 1;
	}

	ถ้า (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		คืน;
	}

#endif
	จำนวนเต็ม buffer_size;
	blasint i, j;

#ifdef SMP
	จำนวนเต็ม nthreads;
#endif


#ifdef SMP
	static จำนวนเต็ม (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, จำนวนเต็ม) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	จำนวนเต็ม (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ถ้า (m == 0)
		คืน;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ถ้า (uplo == 1) {
		สำหรับ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ถ้า (transa & 1) {
				aa = a + lda * i;
			}
			ถ้า (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ถ้า 0
			ถ้า (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			ถ้า (alpha == ZERO)
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

			ถ้า (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			มิฉะนั้น
				nthreads = num_cpu_avail(2);

			ถ้า (nthreads == 1) {
#endif

				ถ้า (!(transa & 1))
				(gemv[(จำนวนเต็ม)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				มิฉะนั้น
				(gemv[(จำนวนเต็ม)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} มิฉะนั้น {
				ถ้า (!(transa & 1))
				(gemv_thread[(จำนวนเต็ม)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				มิฉะนั้น
				(gemv_thread[(จำนวนเต็ม)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} มิฉะนั้น {

		สำหรับ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ถ้า (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ถ้า 0
			ถ้า (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			ถ้า (alpha == ZERO)
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

			ถ้า (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			มิฉะนั้น
				nthreads = num_cpu_avail(2);

			ถ้า (nthreads == 1) {
#endif

				ถ้า (!(transa & 1))
				(gemv[(จำนวนเต็ม)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				มิฉะนั้น
				(gemv[(จำนวนเต็ม)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} มิฉะนั้น {
				ถ้า (!(transa & 1))
				(gemv_thread[(จำนวนเต็ม)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				มิฉะนั้น
				(gemv_thread[(จำนวนเต็ม)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	คืน;
}
