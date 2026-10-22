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

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *UPLO, نویسه *TRANSA, نویسه *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	صحیح transa, transb, uplo;
	blasint info;

	نویسه transA, transB, Uplo;
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

	اگر (transA == 'N')
		transa = 0;
	اگر (transA == 'T')
		transa = 1;

	اگر (transA == 'R')
		transa = 0;
	اگر (transA == 'C')
		transa = 1;

	اگر (transB == 'N')
		transb = 0;
	اگر (transB == 'T')
		transb = 1;

	اگر (transB == 'R')
		transb = 0;
	اگر (transB == 'C')
		transb = 1;

	اگر (Uplo == 'U')
		uplo = 0;
	اگر (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	اگر (transa & 1) nrowa = k;
	nrowb = k;
	اگر (transb & 1) nrowb = m;

	info = 0;

	اگر (ldc < MAX(1, m))
		info = 13;
	اگر (ldb < MAX(1, nrowb))
		info = 10;
	اگر (lda < MAX(1, nrowa))
		info = 8;
	اگر (k < 0)
		info = 5;
	اگر (m < 0)
		info = 4;
	اگر (transb < 0)
		info = 3;
	اگر (transa < 0)
		info = 2;
	اگر (uplo < 0)
		info = 1;

	اگر (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, اندازه(ERROR_NAME));
		بازگشت;
	}
#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order, شمارش CBLAS_UPLO Uplo,
	   شمارش CBLAS_TRANSPOSE TransA, شمارش CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	صحیح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	اگر (order == CblasColMajor) {
		اگر (Uplo == CblasUpper) uplo = 0;
		اگر (Uplo == CblasLower) uplo = 1;

		اگر (TransA == CblasNoTrans)
			transa = 0;
		اگر (TransA == CblasTrans)
			transa = 1;

		اگر (TransA == CblasConjNoTrans)
			transa = 0;
		اگر (TransA == CblasConjTrans)
			transa = 1;

		اگر (TransB == CblasNoTrans)
			transb = 0;
		اگر (TransB == CblasTrans)
			transb = 1;

		اگر (TransB == CblasConjNoTrans)
			transb = 0;
		اگر (TransB == CblasConjTrans)
			transb = 1;

		a = (تهی *)A;
		b = (تهی *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		اگر (transa & 1) nrowa = k;
		nrowb = k;
		اگر (transb & 1)  nrowb = m;

		اگر (ldc < MAX(1, m))
			info = 13;
		اگر (ldb < MAX(1, nrowb))
			info = 10;
		اگر (lda < MAX(1, nrowa))
			info = 8;
		اگر (k < 0)
			info = 5;
		اگر (m < 0)
			info = 4;
		اگر (transb < 0)
			info = 3;
		اگر (transa < 0)
			info = 2;
		اگر (uplo < 0)
			info = 1;
	}

	اگر (order == CblasRowMajor) {

		a = (تهی *)B;
		b = (تهی *)A;

		lda = LDB;
		ldb = LDA;

		اگر (Uplo == CblasUpper) uplo = 0;
		اگر (Uplo == CblasLower) uplo = 1;

		اگر (TransB == CblasNoTrans)
			transa = 0;
		اگر (TransB == CblasTrans)
			transa = 1;

		اگر (TransB == CblasConjNoTrans)
			transa = 0;
		اگر (TransB == CblasConjTrans)
			transa = 1;

		اگر (TransA == CblasNoTrans)
			transb = 0;
		اگر (TransA == CblasTrans)
			transb = 1;

		اگر (TransA == CblasConjNoTrans)
			transb = 0;
		اگر (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		اگر (transa & 1) ncola = k;
		ncolb = k;

		اگر (transb & 1) {
			ncolb = m;
		}

		اگر (ldc < MAX(1,m))
			info = 13;
		اگر (ldb < MAX(1, ncolb))
			info = 8;
		اگر (lda < MAX(1, ncola))
			info = 10;
		اگر (k < 0)
			info = 5;
		اگر (m < 0)
			info = 4;
		اگر (transb < 0)
			info = 2;
		اگر (transa < 0)
			info = 3;
		اگر (uplo < 0)
			info = 1;
	}

	اگر (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, اندازه(ERROR_NAME));
		بازگشت;
	}

#endif
	صحیح buffer_size;
	blasint i, j;

#ifdef SMP
	صحیح nthreads;
#endif


#ifdef SMP
	ایستا صحیح (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, صحیح) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	صحیح (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	اگر (m == 0)
		بازگشت;

	IDEBUG_START;

	ثابت blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	اگر (uplo == 1) {
		برای (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			اگر (transa & 1) {
				aa = a + lda * i;
			}
			اگر (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#اگر 0
			اگر (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, پوچ, 0, پوچ, 0);

			اگر (alpha == ZERO)
				ادامه;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / اندازه(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / اندازه(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			پشته_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			اگر (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			وگرنه
				nthreads = رقم_cpu_avail(2);

			اگر (nthreads == 1) {
#endif

				اگر (!(transa & 1))
				(gemv[(صحیح)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				وگرنه
				(gemv[(صحیح)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} وگرنه {
				اگر (!(transa & 1))
				(gemv_thread[(صحیح)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				وگرنه
				(gemv_thread[(صحیح)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			پشته_FREE(buffer);
		}
	} وگرنه {

		برای (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			اگر (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#اگر 0
			اگر (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, پوچ, 0, پوچ, 0);

			اگر (alpha == ZERO)
				ادامه;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / اندازه(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / اندازه(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			پشته_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			اگر (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			وگرنه
				nthreads = رقم_cpu_avail(2);

			اگر (nthreads == 1) {
#endif

				اگر (!(transa & 1))
				(gemv[(صحیح)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				وگرنه
				(gemv[(صحیح)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} وگرنه {
				اگر (!(transa & 1))
				(gemv_thread[(صحیح)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				وگرنه
				(gemv_thread[(صحیح)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			پشته_FREE(buffer);
		}
	}

	IDEBUG_END;

	بازگشت;
}
