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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, ଅକ୍ଷର *TRANSA, ଅକ୍ଷର *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା transa, transb, uplo;
	blasint info;

	ଅକ୍ଷର transA, transB, Uplo;
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

	ଯଦି (transA == 'N')
		transa = 0;
	ଯଦି (transA == 'T')
		transa = 1;

	ଯଦି (transA == 'R')
		transa = 0;
	ଯଦି (transA == 'C')
		transa = 1;

	ଯଦି (transB == 'N')
		transb = 0;
	ଯଦି (transB == 'T')
		transb = 1;

	ଯଦି (transB == 'R')
		transb = 0;
	ଯଦି (transB == 'C')
		transb = 1;

	ଯଦି (Uplo == 'U')
		uplo = 0;
	ଯଦି (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	ଯଦି (transa & 1) nrowa = k;
	nrowb = k;
	ଯଦି (transb & 1) nrowb = m;

	info = 0;

	ଯଦି (ldc < MAX(1, m))
		info = 13;
	ଯଦି (ldb < MAX(1, nrowb))
		info = 10;
	ଯଦି (lda < MAX(1, nrowa))
		info = 8;
	ଯଦି (k < 0)
		info = 5;
	ଯଦି (m < 0)
		info = 4;
	ଯଦି (transb < 0)
		info = 3;
	ଯଦି (transa < 0)
		info = 2;
	ଯଦି (uplo < 0)
		info = 1;

	ଯଦି (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ଆକାର(ERROR_NAME));
		ଫେରନ୍ତୁ;
	}
#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_UPLO Uplo,
	   ଗଣନା CBLAS_TRANSPOSE TransA, ଗଣନା CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	ପୂର୍ଣ୍ଣସଂଖ୍ୟା transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ଯଦି (order == CblasColMajor) {
		ଯଦି (Uplo == CblasUpper) uplo = 0;
		ଯଦି (Uplo == CblasLower) uplo = 1;

		ଯଦି (TransA == CblasNoTrans)
			transa = 0;
		ଯଦି (TransA == CblasTrans)
			transa = 1;

		ଯଦି (TransA == CblasConjNoTrans)
			transa = 0;
		ଯଦି (TransA == CblasConjTrans)
			transa = 1;

		ଯଦି (TransB == CblasNoTrans)
			transb = 0;
		ଯଦି (TransB == CblasTrans)
			transb = 1;

		ଯଦି (TransB == CblasConjNoTrans)
			transb = 0;
		ଯଦି (TransB == CblasConjTrans)
			transb = 1;

		a = (ଶୂନ୍ୟ *)A;
		b = (ଶୂନ୍ୟ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ଯଦି (transa & 1) nrowa = k;
		nrowb = k;
		ଯଦି (transb & 1)  nrowb = m;

		ଯଦି (ldc < MAX(1, m))
			info = 13;
		ଯଦି (ldb < MAX(1, nrowb))
			info = 10;
		ଯଦି (lda < MAX(1, nrowa))
			info = 8;
		ଯଦି (k < 0)
			info = 5;
		ଯଦି (m < 0)
			info = 4;
		ଯଦି (transb < 0)
			info = 3;
		ଯଦି (transa < 0)
			info = 2;
		ଯଦି (uplo < 0)
			info = 1;
	}

	ଯଦି (order == CblasRowMajor) {

		a = (ଶୂନ୍ୟ *)B;
		b = (ଶୂନ୍ୟ *)A;

		lda = LDB;
		ldb = LDA;

		ଯଦି (Uplo == CblasUpper) uplo = 0;
		ଯଦି (Uplo == CblasLower) uplo = 1;

		ଯଦି (TransB == CblasNoTrans)
			transa = 0;
		ଯଦି (TransB == CblasTrans)
			transa = 1;

		ଯଦି (TransB == CblasConjNoTrans)
			transa = 0;
		ଯଦି (TransB == CblasConjTrans)
			transa = 1;

		ଯଦି (TransA == CblasNoTrans)
			transb = 0;
		ଯଦି (TransA == CblasTrans)
			transb = 1;

		ଯଦି (TransA == CblasConjNoTrans)
			transb = 0;
		ଯଦି (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ଯଦି (transa & 1) ncola = k;
		ncolb = k;

		ଯଦି (transb & 1) {
			ncolb = m;
		}

		ଯଦି (ldc < MAX(1,m))
			info = 13;
		ଯଦି (ldb < MAX(1, ncolb))
			info = 8;
		ଯଦି (lda < MAX(1, ncola))
			info = 10;
		ଯଦି (k < 0)
			info = 5;
		ଯଦି (m < 0)
			info = 4;
		ଯଦି (transb < 0)
			info = 2;
		ଯଦି (transa < 0)
			info = 3;
		ଯଦି (uplo < 0)
			info = 1;
	}

	ଯଦି (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, ଆକାର(ERROR_NAME));
		ଫେରନ୍ତୁ;
	}

#endif
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା buffer_size;
	blasint i, j;

#ifdef SMP
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif


#ifdef SMP
	ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, ପୂର୍ଣ୍ଣସଂଖ୍ୟା) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ଯଦି (m == 0)
		ଫେରନ୍ତୁ;

	IDEBUG_START;

	ଅଚଳ blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ଯଦି (uplo == 1) {
		ପାଇଁ (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ଯଦି (transa & 1) {
				aa = a + lda * i;
			}
			ଯଦି (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ଯଦି 0
			ଯଦି (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ଖାଲି, 0, ଖାଲି, 0);

			ଯଦି (alpha == ZERO)
				ଜାରି;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / ଆକାର(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / ଆକାର(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			ଗଦା_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ଯଦି (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			ଅନ୍ୟଥା
				nthreads = ଅଙ୍କ_cpu_avail(2);

			ଯଦି (nthreads == 1) {
#endif

				ଯଦି (!(transa & 1))
				(gemv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				ଅନ୍ୟଥା
				(gemv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} ଅନ୍ୟଥା {
				ଯଦି (!(transa & 1))
				(gemv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				ଅନ୍ୟଥା
				(gemv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			ଗଦା_FREE(buffer);
		}
	} ଅନ୍ୟଥା {

		ପାଇଁ (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ଯଦି (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ଯଦି 0
			ଯଦି (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ଖାଲି, 0, ଖାଲି, 0);

			ଯଦି (alpha == ZERO)
				ଜାରି;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / ଆକାର(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / ଆକାର(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			ଗଦା_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			ଯଦି (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			ଅନ୍ୟଥା
				nthreads = ଅଙ୍କ_cpu_avail(2);

			ଯଦି (nthreads == 1) {
#endif

				ଯଦି (!(transa & 1))
				(gemv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				ଅନ୍ୟଥା
				(gemv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} ଅନ୍ୟଥା {
				ଯଦି (!(transa & 1))
				(gemv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				ଅନ୍ୟଥା
				(gemv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			ଗଦା_FREE(buffer);
		}
	}

	IDEBUG_END;

	ଫେରନ୍ତୁ;
}
