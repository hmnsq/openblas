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

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANSA, حرف *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	صحيح transa, transb, uplo;
	blasint info;

	حرف transA, transB, Uplo;
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

	إذا (transA == 'N')
		transa = 0;
	إذا (transA == 'T')
		transa = 1;

	إذا (transA == 'R')
		transa = 0;
	إذا (transA == 'C')
		transa = 1;

	إذا (transB == 'N')
		transb = 0;
	إذا (transB == 'T')
		transb = 1;

	إذا (transB == 'R')
		transb = 0;
	إذا (transB == 'C')
		transb = 1;

	إذا (Uplo == 'U')
		uplo = 0;
	إذا (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	إذا (transa & 1) nrowa = k;
	nrowb = k;
	إذا (transb & 1) nrowb = m;

	info = 0;

	إذا (ldc < MAX(1, m))
		info = 13;
	إذا (ldb < MAX(1, nrowb))
		info = 10;
	إذا (lda < MAX(1, nrowa))
		info = 8;
	إذا (k < 0)
		info = 5;
	إذا (m < 0)
		info = 4;
	إذا (transb < 0)
		info = 3;
	إذا (transa < 0)
		info = 2;
	إذا (uplo < 0)
		info = 1;

	إذا (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}
#وإلا

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	صحيح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	إذا (order == CblasColMajor) {
		إذا (Uplo == CblasUpper) uplo = 0;
		إذا (Uplo == CblasLower) uplo = 1;

		إذا (TransA == CblasNoTrans)
			transa = 0;
		إذا (TransA == CblasTrans)
			transa = 1;

		إذا (TransA == CblasConjNoTrans)
			transa = 0;
		إذا (TransA == CblasConjTrans)
			transa = 1;

		إذا (TransB == CblasNoTrans)
			transb = 0;
		إذا (TransB == CblasTrans)
			transb = 1;

		إذا (TransB == CblasConjNoTrans)
			transb = 0;
		إذا (TransB == CblasConjTrans)
			transb = 1;

		a = (فراغ *)A;
		b = (فراغ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		إذا (transa & 1) nrowa = k;
		nrowb = k;
		إذا (transb & 1)  nrowb = m;

		إذا (ldc < MAX(1, m))
			info = 13;
		إذا (ldb < MAX(1, nrowb))
			info = 10;
		إذا (lda < MAX(1, nrowa))
			info = 8;
		إذا (k < 0)
			info = 5;
		إذا (m < 0)
			info = 4;
		إذا (transb < 0)
			info = 3;
		إذا (transa < 0)
			info = 2;
		إذا (uplo < 0)
			info = 1;
	}

	إذا (order == CblasRowMajor) {

		a = (فراغ *)B;
		b = (فراغ *)A;

		lda = LDB;
		ldb = LDA;

		إذا (Uplo == CblasUpper) uplo = 0;
		إذا (Uplo == CblasLower) uplo = 1;

		إذا (TransB == CblasNoTrans)
			transa = 0;
		إذا (TransB == CblasTrans)
			transa = 1;

		إذا (TransB == CblasConjNoTrans)
			transa = 0;
		إذا (TransB == CblasConjTrans)
			transa = 1;

		إذا (TransA == CblasNoTrans)
			transb = 0;
		إذا (TransA == CblasTrans)
			transb = 1;

		إذا (TransA == CblasConjNoTrans)
			transb = 0;
		إذا (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		إذا (transa & 1) ncola = k;
		ncolb = k;

		إذا (transb & 1) {
			ncolb = m;
		}

		إذا (ldc < MAX(1,m))
			info = 13;
		إذا (ldb < MAX(1, ncolb))
			info = 8;
		إذا (lda < MAX(1, ncola))
			info = 10;
		إذا (k < 0)
			info = 5;
		إذا (m < 0)
			info = 4;
		إذا (transb < 0)
			info = 2;
		إذا (transa < 0)
			info = 3;
		إذا (uplo < 0)
			info = 1;
	}

	إذا (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}

#endif
	صحيح buffer_size;
	blasint i, j;

#ifdef SMP
	صحيح nthreads;
#endif


#ifdef SMP
	ساكن صحيح (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, صحيح) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	صحيح (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	إذا (m == 0)
		أرجع;

	IDEBUG_START;

	ثابت blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	إذا (uplo == 1) {
		لكل (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			إذا (transa & 1) {
				aa = a + lda * i;
			}
			إذا (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#إذا 0
			إذا (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, عدم, 0, عدم, 0);

			إذا (alpha == ZERO)
				استمر;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / حجم(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / حجم(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			کەڵەکە_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			إذا (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			وإلا
				nthreads = ژمارە_cpu_avail(2);

			إذا (nthreads == 1) {
#endif

				إذا (!(transa & 1))
				(gemv[(صحيح)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				وإلا
				(gemv[(صحيح)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} وإلا {
				إذا (!(transa & 1))
				(gemv_thread[(صحيح)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				وإلا
				(gemv_thread[(صحيح)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			کەڵەکە_FREE(buffer);
		}
	} وإلا {

		لكل (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			إذا (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#إذا 0
			إذا (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, عدم, 0, عدم, 0);

			إذا (alpha == ZERO)
				استمر;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / حجم(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / حجم(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			کەڵەکە_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			إذا (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			وإلا
				nthreads = ژمارە_cpu_avail(2);

			إذا (nthreads == 1) {
#endif

				إذا (!(transa & 1))
				(gemv[(صحيح)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				وإلا
				(gemv[(صحيح)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} وإلا {
				إذا (!(transa & 1))
				(gemv_thread[(صحيح)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				وإلا
				(gemv_thread[(صحيح)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			کەڵەکە_FREE(buffer);
		}
	}

	IDEBUG_END;

	أرجع;
}
