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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <stdlib.h>
#அந்தர்பாவ "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, எழுத்து *TRANSA, எழுத்து *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	முழு transa, transb, uplo;
	blasint info;

	எழுத்து transA, transB, Uplo;
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

	எனில் (transA == 'N')
		transa = 0;
	எனில் (transA == 'T')
		transa = 1;

	எனில் (transA == 'R')
		transa = 0;
	எனில் (transA == 'C')
		transa = 1;

	எனில் (transB == 'N')
		transb = 0;
	எனில் (transB == 'T')
		transb = 1;

	எனில் (transB == 'R')
		transb = 0;
	எனில் (transB == 'C')
		transb = 1;

	எனில் (Uplo == 'U')
		uplo = 0;
	எனில் (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	எனில் (transa & 1) nrowa = k;
	nrowb = k;
	எனில் (transb & 1) nrowb = m;

	info = 0;

	எனில் (ldc < MAX(1, m))
		info = 13;
	எனில் (ldb < MAX(1, nrowb))
		info = 10;
	எனில் (lda < MAX(1, nrowa))
		info = 8;
	எனில் (k < 0)
		info = 5;
	எனில் (m < 0)
		info = 4;
	எனில் (transb < 0)
		info = 3;
	எனில் (transa < 0)
		info = 2;
	எனில் (uplo < 0)
		info = 1;

	எனில் (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, அளவு(ERROR_NAME));
		திரும்பு;
	}
#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_UPLO Uplo,
	   எண்ணிக்கை CBLAS_TRANSPOSE TransA, எண்ணிக்கை CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	முழு transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	எனில் (order == CblasColMajor) {
		எனில் (Uplo == CblasUpper) uplo = 0;
		எனில் (Uplo == CblasLower) uplo = 1;

		எனில் (TransA == CblasNoTrans)
			transa = 0;
		எனில் (TransA == CblasTrans)
			transa = 1;

		எனில் (TransA == CblasConjNoTrans)
			transa = 0;
		எனில் (TransA == CblasConjTrans)
			transa = 1;

		எனில் (TransB == CblasNoTrans)
			transb = 0;
		எனில் (TransB == CblasTrans)
			transb = 1;

		எனில் (TransB == CblasConjNoTrans)
			transb = 0;
		எனில் (TransB == CblasConjTrans)
			transb = 1;

		a = (வெறுமை *)A;
		b = (வெறுமை *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		எனில் (transa & 1) nrowa = k;
		nrowb = k;
		எனில் (transb & 1)  nrowb = m;

		எனில் (ldc < MAX(1, m))
			info = 13;
		எனில் (ldb < MAX(1, nrowb))
			info = 10;
		எனில் (lda < MAX(1, nrowa))
			info = 8;
		எனில் (k < 0)
			info = 5;
		எனில் (m < 0)
			info = 4;
		எனில் (transb < 0)
			info = 3;
		எனில் (transa < 0)
			info = 2;
		எனில் (uplo < 0)
			info = 1;
	}

	எனில் (order == CblasRowMajor) {

		a = (வெறுமை *)B;
		b = (வெறுமை *)A;

		lda = LDB;
		ldb = LDA;

		எனில் (Uplo == CblasUpper) uplo = 0;
		எனில் (Uplo == CblasLower) uplo = 1;

		எனில் (TransB == CblasNoTrans)
			transa = 0;
		எனில் (TransB == CblasTrans)
			transa = 1;

		எனில் (TransB == CblasConjNoTrans)
			transa = 0;
		எனில் (TransB == CblasConjTrans)
			transa = 1;

		எனில் (TransA == CblasNoTrans)
			transb = 0;
		எனில் (TransA == CblasTrans)
			transb = 1;

		எனில் (TransA == CblasConjNoTrans)
			transb = 0;
		எனில் (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		எனில் (transa & 1) ncola = k;
		ncolb = k;

		எனில் (transb & 1) {
			ncolb = m;
		}

		எனில் (ldc < MAX(1,m))
			info = 13;
		எனில் (ldb < MAX(1, ncolb))
			info = 8;
		எனில் (lda < MAX(1, ncola))
			info = 10;
		எனில் (k < 0)
			info = 5;
		எனில் (m < 0)
			info = 4;
		எனில் (transb < 0)
			info = 2;
		எனில் (transa < 0)
			info = 3;
		எனில் (uplo < 0)
			info = 1;
	}

	எனில் (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, அளவு(ERROR_NAME));
		திரும்பு;
	}

#endif
	முழு buffer_size;
	blasint i, j;

#ifdef SMP
	முழு nthreads;
#endif


#ifdef SMP
	நிலையான முழு (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, முழு) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	முழு (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	எனில் (m == 0)
		திரும்பு;

	IDEBUG_START;

	மாறிலி blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	எனில் (uplo == 1) {
		க்கு (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			எனில் (transa & 1) {
				aa = a + lda * i;
			}
			எனில் (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#எனில் 0
			எனில் (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, வெற்று, 0, வெற்று, 0);

			எனில் (alpha == ZERO)
				தொடர்;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / அளவு(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / அளவு(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			அடுக்கு_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			எனில் (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			இல்லையெனில்
				nthreads = இலக்கம்_cpu_avail(2);

			எனில் (nthreads == 1) {
#endif

				எனில் (!(transa & 1))
				(gemv[(முழு)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				இல்லையெனில்
				(gemv[(முழு)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} இல்லையெனில் {
				எனில் (!(transa & 1))
				(gemv_thread[(முழு)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				இல்லையெனில்
				(gemv_thread[(முழு)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			அடுக்கு_FREE(buffer);
		}
	} இல்லையெனில் {

		க்கு (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			எனில் (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#எனில் 0
			எனில் (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, வெற்று, 0, வெற்று, 0);

			எனில் (alpha == ZERO)
				தொடர்;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / அளவு(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / அளவு(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			அடுக்கு_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			எனில் (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			இல்லையெனில்
				nthreads = இலக்கம்_cpu_avail(2);

			எனில் (nthreads == 1) {
#endif

				எனில் (!(transa & 1))
				(gemv[(முழு)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				இல்லையெனில்
				(gemv[(முழு)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} இல்லையெனில் {
				எனில் (!(transa & 1))
				(gemv_thread[(முழு)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				இல்லையெனில்
				(gemv_thread[(முழு)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			அடுக்கு_FREE(buffer);
		}
	}

	IDEBUG_END;

	திரும்பு;
}
