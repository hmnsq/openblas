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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <stdlib.h>
#অন্তর্ভাব "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(অক্ষর *UPLO, অক্ষর *TRANSA, অক্ষর *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	পূর্ণসংখ্যা transa, transb, uplo;
	blasint info;

	অক্ষর transA, transB, Uplo;
	blasint nrowa, nrowb;
	IFLOAT *বাফার;
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

	যদি (transA == 'N')
		transa = 0;
	যদি (transA == 'T')
		transa = 1;

	যদি (transA == 'R')
		transa = 0;
	যদি (transA == 'C')
		transa = 1;

	যদি (transB == 'N')
		transb = 0;
	যদি (transB == 'T')
		transb = 1;

	যদি (transB == 'R')
		transb = 0;
	যদি (transB == 'C')
		transb = 1;

	যদি (Uplo == 'U')
		uplo = 0;
	যদি (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	যদি (transa & 1) nrowa = k;
	nrowb = k;
	যদি (transb & 1) nrowb = m;

	info = 0;

	যদি (ldc < MAX(1, m))
		info = 13;
	যদি (ldb < MAX(1, nrowb))
		info = 10;
	যদি (lda < MAX(1, nrowa))
		info = 8;
	যদি (k < 0)
		info = 5;
	যদি (m < 0)
		info = 4;
	যদি (transb < 0)
		info = 3;
	যদি (transa < 0)
		info = 2;
	যদি (uplo < 0)
		info = 1;

	যদি (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, আকার(ERROR_NAME));
		ফেরত;
	}
#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo,
	   গণনা CBLAS_TRANSPOSE TransA, গণনা CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	পূর্ণসংখ্যা transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *বাফার;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	যদি (order == CblasColMajor) {
		যদি (Uplo == CblasUpper) uplo = 0;
		যদি (Uplo == CblasLower) uplo = 1;

		যদি (TransA == CblasNoTrans)
			transa = 0;
		যদি (TransA == CblasTrans)
			transa = 1;

		যদি (TransA == CblasConjNoTrans)
			transa = 0;
		যদি (TransA == CblasConjTrans)
			transa = 1;

		যদি (TransB == CblasNoTrans)
			transb = 0;
		যদি (TransB == CblasTrans)
			transb = 1;

		যদি (TransB == CblasConjNoTrans)
			transb = 0;
		যদি (TransB == CblasConjTrans)
			transb = 1;

		a = (শূন্য *)A;
		b = (শূন্য *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		যদি (transa & 1) nrowa = k;
		nrowb = k;
		যদি (transb & 1)  nrowb = m;

		যদি (ldc < MAX(1, m))
			info = 13;
		যদি (ldb < MAX(1, nrowb))
			info = 10;
		যদি (lda < MAX(1, nrowa))
			info = 8;
		যদি (k < 0)
			info = 5;
		যদি (m < 0)
			info = 4;
		যদি (transb < 0)
			info = 3;
		যদি (transa < 0)
			info = 2;
		যদি (uplo < 0)
			info = 1;
	}

	যদি (order == CblasRowMajor) {

		a = (শূন্য *)B;
		b = (শূন্য *)A;

		lda = LDB;
		ldb = LDA;

		যদি (Uplo == CblasUpper) uplo = 0;
		যদি (Uplo == CblasLower) uplo = 1;

		যদি (TransB == CblasNoTrans)
			transa = 0;
		যদি (TransB == CblasTrans)
			transa = 1;

		যদি (TransB == CblasConjNoTrans)
			transa = 0;
		যদি (TransB == CblasConjTrans)
			transa = 1;

		যদি (TransA == CblasNoTrans)
			transb = 0;
		যদি (TransA == CblasTrans)
			transb = 1;

		যদি (TransA == CblasConjNoTrans)
			transb = 0;
		যদি (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		যদি (transa & 1) ncola = k;
		ncolb = k;

		যদি (transb & 1) {
			ncolb = m;
		}

		যদি (ldc < MAX(1,m))
			info = 13;
		যদি (ldb < MAX(1, ncolb))
			info = 8;
		যদি (lda < MAX(1, ncola))
			info = 10;
		যদি (k < 0)
			info = 5;
		যদি (m < 0)
			info = 4;
		যদি (transb < 0)
			info = 2;
		যদি (transa < 0)
			info = 3;
		যদি (uplo < 0)
			info = 1;
	}

	যদি (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, আকার(ERROR_NAME));
		ফেরত;
	}

#endif
	পূর্ণসংখ্যা বাফার_আকার;
	blasint i, j;

#ifdef SMP
	পূর্ণসংখ্যা nthreads;
#endif


#ifdef SMP
	স্থির পূর্ণসংখ্যা (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, পূর্ণসংখ্যা) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	পূর্ণসংখ্যা (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	যদি (m == 0)
		ফেরত;

	IDEBUG_START;

	ধ্রুবক blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	যদি (uplo == 1) {
		জন্য (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			যদি (transa & 1) {
				aa = a + lda * i;
			}
			যদি (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#যদি 0
			যদি (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ফাঁকা, 0, ফাঁকা, 0);

			যদি (alpha == ZERO)
				চালু;
#endif

			IDEBUG_START;

			বাফার_আকার = j + k + 128 / আকার(FLOAT);
#ifdef WINDOWS_ABI
			বাফার_আকার += 160 / আকার(FLOAT);
#endif
			// for alignment
			বাফার_আকার = (বাফার_আকার + 3) & ~3;
			স্তূপ_ALLOC(বাফার_আকার, IFLOAT, বাফার);

#ifdef SMP

			যদি (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			অন্যথা
				nthreads = অঙ্ক_cpu_avail(2);

			যদি (nthreads == 1) {
#endif

				যদি (!(transa & 1))
				(gemv[(পূর্ণসংখ্যা)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				অন্যথা
				(gemv[(পূর্ণসংখ্যা)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} অন্যথা {
				যদি (!(transa & 1))
				(gemv_thread[(পূর্ণসংখ্যা)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				অন্যথা
				(gemv_thread[(পূর্ণসংখ্যা)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			স্তূপ_FREE(বাফার);
		}
	} অন্যথা {

		জন্য (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			যদি (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#যদি 0
			যদি (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ফাঁকা, 0, ফাঁকা, 0);

			যদি (alpha == ZERO)
				চালু;
#endif
			IDEBUG_START;

			বাফার_আকার = j + k + 128 / আকার(FLOAT);
#ifdef WINDOWS_ABI
			বাফার_আকার += 160 / আকার(FLOAT);
#endif
			// for alignment
			বাফার_আকার = (বাফার_আকার + 3) & ~3;
			স্তূপ_ALLOC(বাফার_আকার, IFLOAT, বাফার);

#ifdef SMP

			যদি (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			অন্যথা
				nthreads = অঙ্ক_cpu_avail(2);

			যদি (nthreads == 1) {
#endif

				যদি (!(transa & 1))
				(gemv[(পূর্ণসংখ্যা)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				অন্যথা
				(gemv[(পূর্ণসংখ্যা)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} অন্যথা {
				যদি (!(transa & 1))
				(gemv_thread[(পূর্ণসংখ্যা)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				অন্যথা
				(gemv_thread[(পূর্ণসংখ্যা)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			স্তূপ_FREE(বাফার);
		}
	}

	IDEBUG_END;

	ফেরত;
}
