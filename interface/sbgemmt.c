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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANSA, वर्ण *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	पूर्णाङ्क transa, transb, uplo;
	blasint info;

	वर्ण transA, transB, Uplo;
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

	जइ (transA == 'N')
		transa = 0;
	जइ (transA == 'T')
		transa = 1;

	जइ (transA == 'R')
		transa = 0;
	जइ (transA == 'C')
		transa = 1;

	जइ (transB == 'N')
		transb = 0;
	जइ (transB == 'T')
		transb = 1;

	जइ (transB == 'R')
		transb = 0;
	जइ (transB == 'C')
		transb = 1;

	जइ (Uplo == 'U')
		uplo = 0;
	जइ (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	जइ (transa & 1) nrowa = k;
	nrowb = k;
	जइ (transb & 1) nrowb = m;

	info = 0;

	जइ (ldc < MAX(1, m))
		info = 13;
	जइ (ldb < MAX(1, nrowb))
		info = 10;
	जइ (lda < MAX(1, nrowa))
		info = 8;
	जइ (k < 0)
		info = 5;
	जइ (m < 0)
		info = 4;
	जइ (transb < 0)
		info = 3;
	जइ (transa < 0)
		info = 2;
	जइ (uplo < 0)
		info = 1;

	जइ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, परिमाण(ERROR_NAME));
		प्रत्यागम;
	}
#अण्णहा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	पूर्णाङ्क transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	जइ (order == CblasColMajor) {
		जइ (Uplo == CblasUpper) uplo = 0;
		जइ (Uplo == CblasLower) uplo = 1;

		जइ (TransA == CblasNoTrans)
			transa = 0;
		जइ (TransA == CblasTrans)
			transa = 1;

		जइ (TransA == CblasConjNoTrans)
			transa = 0;
		जइ (TransA == CblasConjTrans)
			transa = 1;

		जइ (TransB == CblasNoTrans)
			transb = 0;
		जइ (TransB == CblasTrans)
			transb = 1;

		जइ (TransB == CblasConjNoTrans)
			transb = 0;
		जइ (TransB == CblasConjTrans)
			transb = 1;

		a = (शून्य *)A;
		b = (शून्य *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		जइ (transa & 1) nrowa = k;
		nrowb = k;
		जइ (transb & 1)  nrowb = m;

		जइ (ldc < MAX(1, m))
			info = 13;
		जइ (ldb < MAX(1, nrowb))
			info = 10;
		जइ (lda < MAX(1, nrowa))
			info = 8;
		जइ (k < 0)
			info = 5;
		जइ (m < 0)
			info = 4;
		जइ (transb < 0)
			info = 3;
		जइ (transa < 0)
			info = 2;
		जइ (uplo < 0)
			info = 1;
	}

	जइ (order == CblasRowMajor) {

		a = (शून्य *)B;
		b = (शून्य *)A;

		lda = LDB;
		ldb = LDA;

		जइ (Uplo == CblasUpper) uplo = 0;
		जइ (Uplo == CblasLower) uplo = 1;

		जइ (TransB == CblasNoTrans)
			transa = 0;
		जइ (TransB == CblasTrans)
			transa = 1;

		जइ (TransB == CblasConjNoTrans)
			transa = 0;
		जइ (TransB == CblasConjTrans)
			transa = 1;

		जइ (TransA == CblasNoTrans)
			transb = 0;
		जइ (TransA == CblasTrans)
			transb = 1;

		जइ (TransA == CblasConjNoTrans)
			transb = 0;
		जइ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		जइ (transa & 1) ncola = k;
		ncolb = k;

		जइ (transb & 1) {
			ncolb = m;
		}

		जइ (ldc < MAX(1,m))
			info = 13;
		जइ (ldb < MAX(1, ncolb))
			info = 8;
		जइ (lda < MAX(1, ncola))
			info = 10;
		जइ (k < 0)
			info = 5;
		जइ (m < 0)
			info = 4;
		जइ (transb < 0)
			info = 2;
		जइ (transa < 0)
			info = 3;
		जइ (uplo < 0)
			info = 1;
	}

	जइ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, परिमाण(ERROR_NAME));
		प्रत्यागम;
	}

#endif
	पूर्णाङ्क buffer_size;
	blasint i, j;

#ifdef SMP
	पूर्णाङ्क nthreads;
#endif


#ifdef SMP
	स्थिर पूर्णाङ्क (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, पूर्णाङ्क) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	पूर्णाङ्क (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	जइ (m == 0)
		प्रत्यागम;

	IDEBUG_START;

	ध्रुव blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	जइ (uplo == 1) {
		प्रति (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			जइ (transa & 1) {
				aa = a + lda * i;
			}
			जइ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#जइ 0
			जइ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, रिक्त, 0, रिक्त, 0);

			जइ (alpha == ZERO)
				अनुवर्त;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / परिमाण(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / परिमाण(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			राशिः_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			जइ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			अण्णहा
				nthreads = अङ्कः_cpu_avail(2);

			जइ (nthreads == 1) {
#endif

				जइ (!(transa & 1))
				(gemv[(पूर्णाङ्क)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				अण्णहा
				(gemv[(पूर्णाङ्क)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} अण्णहा {
				जइ (!(transa & 1))
				(gemv_thread[(पूर्णाङ्क)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				अण्णहा
				(gemv_thread[(पूर्णाङ्क)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			राशिः_FREE(buffer);
		}
	} अण्णहा {

		प्रति (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			जइ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#जइ 0
			जइ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, रिक्त, 0, रिक्त, 0);

			जइ (alpha == ZERO)
				अनुवर्त;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / परिमाण(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / परिमाण(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			राशिः_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			जइ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			अण्णहा
				nthreads = अङ्कः_cpu_avail(2);

			जइ (nthreads == 1) {
#endif

				जइ (!(transa & 1))
				(gemv[(पूर्णाङ्क)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				अण्णहा
				(gemv[(पूर्णाङ्क)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} अण्णहा {
				जइ (!(transa & 1))
				(gemv_thread[(पूर्णाङ्क)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				अण्णहा
				(gemv_thread[(पूर्णाङ्क)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			राशिः_FREE(buffer);
		}
	}

	IDEBUG_END;

	प्रत्यागम;
}
