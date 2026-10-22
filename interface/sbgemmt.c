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

	यदि (transA == 'N')
		transa = 0;
	यदि (transA == 'T')
		transa = 1;

	यदि (transA == 'R')
		transa = 0;
	यदि (transA == 'C')
		transa = 1;

	यदि (transB == 'N')
		transb = 0;
	यदि (transB == 'T')
		transb = 1;

	यदि (transB == 'R')
		transb = 0;
	यदि (transB == 'C')
		transb = 1;

	यदि (Uplo == 'U')
		uplo = 0;
	यदि (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	यदि (transa & 1) nrowa = k;
	nrowb = k;
	यदि (transb & 1) nrowb = m;

	info = 0;

	यदि (ldc < MAX(1, m))
		info = 13;
	यदि (ldb < MAX(1, nrowb))
		info = 10;
	यदि (lda < MAX(1, nrowa))
		info = 8;
	यदि (k < 0)
		info = 5;
	यदि (m < 0)
		info = 4;
	यदि (transb < 0)
		info = 3;
	यदि (transa < 0)
		info = 2;
	यदि (uplo < 0)
		info = 1;

	यदि (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, परिमाण(ERROR_NAME));
		प्रत्यागम;
	}
#अन्यथा

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

	यदि (order == CblasColMajor) {
		यदि (Uplo == CblasUpper) uplo = 0;
		यदि (Uplo == CblasLower) uplo = 1;

		यदि (TransA == CblasNoTrans)
			transa = 0;
		यदि (TransA == CblasTrans)
			transa = 1;

		यदि (TransA == CblasConjNoTrans)
			transa = 0;
		यदि (TransA == CblasConjTrans)
			transa = 1;

		यदि (TransB == CblasNoTrans)
			transb = 0;
		यदि (TransB == CblasTrans)
			transb = 1;

		यदि (TransB == CblasConjNoTrans)
			transb = 0;
		यदि (TransB == CblasConjTrans)
			transb = 1;

		a = (शून्य *)A;
		b = (शून्य *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		यदि (transa & 1) nrowa = k;
		nrowb = k;
		यदि (transb & 1)  nrowb = m;

		यदि (ldc < MAX(1, m))
			info = 13;
		यदि (ldb < MAX(1, nrowb))
			info = 10;
		यदि (lda < MAX(1, nrowa))
			info = 8;
		यदि (k < 0)
			info = 5;
		यदि (m < 0)
			info = 4;
		यदि (transb < 0)
			info = 3;
		यदि (transa < 0)
			info = 2;
		यदि (uplo < 0)
			info = 1;
	}

	यदि (order == CblasRowMajor) {

		a = (शून्य *)B;
		b = (शून्य *)A;

		lda = LDB;
		ldb = LDA;

		यदि (Uplo == CblasUpper) uplo = 0;
		यदि (Uplo == CblasLower) uplo = 1;

		यदि (TransB == CblasNoTrans)
			transa = 0;
		यदि (TransB == CblasTrans)
			transa = 1;

		यदि (TransB == CblasConjNoTrans)
			transa = 0;
		यदि (TransB == CblasConjTrans)
			transa = 1;

		यदि (TransA == CblasNoTrans)
			transb = 0;
		यदि (TransA == CblasTrans)
			transb = 1;

		यदि (TransA == CblasConjNoTrans)
			transb = 0;
		यदि (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		यदि (transa & 1) ncola = k;
		ncolb = k;

		यदि (transb & 1) {
			ncolb = m;
		}

		यदि (ldc < MAX(1,m))
			info = 13;
		यदि (ldb < MAX(1, ncolb))
			info = 8;
		यदि (lda < MAX(1, ncola))
			info = 10;
		यदि (k < 0)
			info = 5;
		यदि (m < 0)
			info = 4;
		यदि (transb < 0)
			info = 2;
		यदि (transa < 0)
			info = 3;
		यदि (uplo < 0)
			info = 1;
	}

	यदि (info >= 0) {
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


	यदि (m == 0)
		प्रत्यागम;

	IDEBUG_START;

	ध्रुव blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	यदि (uplo == 1) {
		प्रति (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			यदि (transa & 1) {
				aa = a + lda * i;
			}
			यदि (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#यदि 0
			यदि (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, रिक्त, 0, रिक्त, 0);

			यदि (alpha == ZERO)
				अनुवर्त;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / परिमाण(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / परिमाण(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			चळत_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			यदि (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			अन्यथा
				nthreads = अंक_cpu_avail(2);

			यदि (nthreads == 1) {
#endif

				यदि (!(transa & 1))
				(gemv[(पूर्णाङ्क)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				अन्यथा
				(gemv[(पूर्णाङ्क)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} अन्यथा {
				यदि (!(transa & 1))
				(gemv_thread[(पूर्णाङ्क)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				अन्यथा
				(gemv_thread[(पूर्णाङ्क)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			चळत_FREE(buffer);
		}
	} अन्यथा {

		प्रति (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			यदि (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#यदि 0
			यदि (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, रिक्त, 0, रिक्त, 0);

			यदि (alpha == ZERO)
				अनुवर्त;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / परिमाण(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / परिमाण(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			चळत_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			यदि (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			अन्यथा
				nthreads = अंक_cpu_avail(2);

			यदि (nthreads == 1) {
#endif

				यदि (!(transa & 1))
				(gemv[(पूर्णाङ्क)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				अन्यथा
				(gemv[(पूर्णाङ्क)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} अन्यथा {
				यदि (!(transa & 1))
				(gemv_thread[(पूर्णाङ्क)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				अन्यथा
				(gemv_thread[(पूर्णाङ्क)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			चळत_FREE(buffer);
		}
	}

	IDEBUG_END;

	प्रत्यागम;
}
