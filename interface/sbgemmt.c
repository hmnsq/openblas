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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ <stdlib.h>
#అన్తర్భావ "common.h"

#నిర్వచించు SMP_పరిమితి_MIN 65536.0
#నిర్వచించు ERROR_NAME "SBGEMMT "

#నిర్వచితం_కాకపోతే GEMM_MULTITHREAD_పరిమితి
#నిర్వచించు GEMM_MULTITHREAD_పరిమితి 4
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, అక్షరం *TRANSA, అక్షరం *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	పూర్ణాంకం transa, transb, uplo;
	blasint info;

	అక్షరం transA, transB, Uplo;
	blasint nrowa, nrowb;
	IFLOAT *బఫరు;
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

	ఒకవేళ (transA == 'N')
		transa = 0;
	ఒకవేళ (transA == 'T')
		transa = 1;

	ఒకవేళ (transA == 'R')
		transa = 0;
	ఒకవేళ (transA == 'C')
		transa = 1;

	ఒకవేళ (transB == 'N')
		transb = 0;
	ఒకవేళ (transB == 'T')
		transb = 1;

	ఒకవేళ (transB == 'R')
		transb = 0;
	ఒకవేళ (transB == 'C')
		transb = 1;

	ఒకవేళ (Uplo == 'U')
		uplo = 0;
	ఒకవేళ (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	ఒకవేళ (transa & 1) nrowa = k;
	nrowb = k;
	ఒకవేళ (transb & 1) nrowb = m;

	info = 0;

	ఒకవేళ (ldc < MAX(1, m))
		info = 13;
	ఒకవేళ (ldb < MAX(1, nrowb))
		info = 10;
	ఒకవేళ (lda < MAX(1, nrowa))
		info = 8;
	ఒకవేళ (k < 0)
		info = 5;
	ఒకవేళ (m < 0)
		info = 4;
	ఒకవేళ (transb < 0)
		info = 3;
	ఒకవేళ (transa < 0)
		info = 2;
	ఒకవేళ (uplo < 0)
		info = 1;

	ఒకవేళ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
		తిరిగి;
	}
#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_UPLO Uplo,
	   గణన CBLAS_TRANSPOSE TransA, గణన CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	పూర్ణాంకం transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *బఫరు;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
		ఒకవేళ (Uplo == Cblas_పై) uplo = 0;
		ఒకవేళ (Uplo == Cblas_క్రింది) uplo = 1;

		ఒకవేళ (TransA == CblasNoTrans)
			transa = 0;
		ఒకవేళ (TransA == CblasTrans)
			transa = 1;

		ఒకవేళ (TransA == CblasConjNoTrans)
			transa = 0;
		ఒకవేళ (TransA == CblasConjTrans)
			transa = 1;

		ఒకవేళ (TransB == CblasNoTrans)
			transb = 0;
		ఒకవేళ (TransB == CblasTrans)
			transb = 1;

		ఒకవేళ (TransB == CblasConjNoTrans)
			transb = 0;
		ఒకవేళ (TransB == CblasConjTrans)
			transb = 1;

		a = (శూన్య *)A;
		b = (శూన్య *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		ఒకవేళ (transa & 1) nrowa = k;
		nrowb = k;
		ఒకవేళ (transb & 1)  nrowb = m;

		ఒకవేళ (ldc < MAX(1, m))
			info = 13;
		ఒకవేళ (ldb < MAX(1, nrowb))
			info = 10;
		ఒకవేళ (lda < MAX(1, nrowa))
			info = 8;
		ఒకవేళ (k < 0)
			info = 5;
		ఒకవేళ (m < 0)
			info = 4;
		ఒకవేళ (transb < 0)
			info = 3;
		ఒకవేళ (transa < 0)
			info = 2;
		ఒకవేళ (uplo < 0)
			info = 1;
	}

	ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {

		a = (శూన్య *)B;
		b = (శూన్య *)A;

		lda = LDB;
		ldb = LDA;

		ఒకవేళ (Uplo == Cblas_పై) uplo = 0;
		ఒకవేళ (Uplo == Cblas_క్రింది) uplo = 1;

		ఒకవేళ (TransB == CblasNoTrans)
			transa = 0;
		ఒకవేళ (TransB == CblasTrans)
			transa = 1;

		ఒకవేళ (TransB == CblasConjNoTrans)
			transa = 0;
		ఒకవేళ (TransB == CblasConjTrans)
			transa = 1;

		ఒకవేళ (TransA == CblasNoTrans)
			transb = 0;
		ఒకవేళ (TransA == CblasTrans)
			transb = 1;

		ఒకవేళ (TransA == CblasConjNoTrans)
			transb = 0;
		ఒకవేళ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		ఒకవేళ (transa & 1) ncola = k;
		ncolb = k;

		ఒకవేళ (transb & 1) {
			ncolb = m;
		}

		ఒకవేళ (ldc < MAX(1,m))
			info = 13;
		ఒకవేళ (ldb < MAX(1, ncolb))
			info = 8;
		ఒకవేళ (lda < MAX(1, ncola))
			info = 10;
		ఒకవేళ (k < 0)
			info = 5;
		ఒకవేళ (m < 0)
			info = 4;
		ఒకవేళ (transb < 0)
			info = 2;
		ఒకవేళ (transa < 0)
			info = 3;
		ఒకవేళ (uplo < 0)
			info = 1;
	}

	ఒకవేళ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
		తిరిగి;
	}

#అయితే_ముగింపు
	పూర్ణాంకం బఫరు_కొలత;
	blasint i, j;

#నిర్వచితమైతే SMP
	పూర్ణాంకం nthreads;
#అయితే_ముగింపు


#నిర్వచితమైతే SMP
	స్థిర పూర్ణాంకం (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, పూర్ణాంకం) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#అయితే_ముగింపు
	పూర్ణాంకం (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	ఒకవేళ (m == 0)
		తిరిగి;

	IDEBUG_ప్రారంభం;

	స్థిరరాశి blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	ఒకవేళ (uplo == 1) {
		కొరకు (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			ఒకవేళ (transa & 1) {
				aa = a + lda * i;
			}
			ఒకవేళ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#ఒకవేళ 0
			ఒకవేళ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ఖాళీ, 0, ఖాళీ, 0);

			ఒకవేళ (alpha == ZERO)
				కొనసాగు;
#అయితే_ముగింపు

			IDEBUG_ప్రారంభం;

			బఫరు_కొలత = j + k + 128 / పరిమాణం(FLOAT);
#నిర్వచితమైతే WINDOWS_ABI
			బఫరు_కొలత += 160 / పరిమాణం(FLOAT);
#అయితే_ముగింపు
			// for alignment
			బఫరు_కొలత = (బఫరు_కొలత + 3) & ~3;
			దొంతర_ALLOC(బఫరు_కొలత, IFLOAT, బఫరు);

#నిర్వచితమైతే SMP

			ఒకవేళ (1L * j * k < 2304L * GEMM_MULTITHREAD_పరిమితి)
				nthreads = 1;
			లేకపోతే
				nthreads = అంకె_cpu_avail(2);

			ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

				ఒకవేళ (!(transa & 1))
				(gemv[(పూర్ణాంకం)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				లేకపోతే
				(gemv[(పూర్ణాంకం)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#నిర్వచితమైతే SMP
			} లేకపోతే {
				ఒకవేళ (!(transa & 1))
				(gemv_thread[(పూర్ణాంకం)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				లేకపోతే
				(gemv_thread[(పూర్ణాంకం)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#అయితే_ముగింపు

			దొంతర_FREE(బఫరు);
		}
	} లేకపోతే {

		కొరకు (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			ఒకవేళ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#ఒకవేళ 0
			ఒకవేళ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ఖాళీ, 0, ఖాళీ, 0);

			ఒకవేళ (alpha == ZERO)
				కొనసాగు;
#అయితే_ముగింపు
			IDEBUG_ప్రారంభం;

			బఫరు_కొలత = j + k + 128 / పరిమాణం(FLOAT);
#నిర్వచితమైతే WINDOWS_ABI
			బఫరు_కొలత += 160 / పరిమాణం(FLOAT);
#అయితే_ముగింపు
			// for alignment
			బఫరు_కొలత = (బఫరు_కొలత + 3) & ~3;
			దొంతర_ALLOC(బఫరు_కొలత, IFLOAT, బఫరు);

#నిర్వచితమైతే SMP

			ఒకవేళ (1L * j * k < 2304L * GEMM_MULTITHREAD_పరిమితి)
				nthreads = 1;
			లేకపోతే
				nthreads = అంకె_cpu_avail(2);

			ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

				ఒకవేళ (!(transa & 1))
				(gemv[(పూర్ణాంకం)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				లేకపోతే
				(gemv[(పూర్ణాంకం)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#నిర్వచితమైతే SMP
			} లేకపోతే {
				ఒకవేళ (!(transa & 1))
				(gemv_thread[(పూర్ణాంకం)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				లేకపోతే
				(gemv_thread[(పూర్ణాంకం)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#అయితే_ముగింపు

			దొంతర_FREE(బఫరు);
		}
	}

	IDEBUG_చివర;

	తిరిగి;
}
