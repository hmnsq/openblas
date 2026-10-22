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
	පූර්ණ_සංඛ්‍යා transa, transb, uplo;
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

	නම් (transA == 'N')
		transa = 0;
	නම් (transA == 'T')
		transa = 1;

	නම් (transA == 'R')
		transa = 0;
	නම් (transA == 'C')
		transa = 1;

	නම් (transB == 'N')
		transb = 0;
	නම් (transB == 'T')
		transb = 1;

	නම් (transB == 'R')
		transb = 0;
	නම් (transB == 'C')
		transb = 1;

	නම් (Uplo == 'U')
		uplo = 0;
	නම් (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	නම් (transa & 1) nrowa = k;
	nrowb = k;
	නම් (transb & 1) nrowb = m;

	info = 0;

	නම් (ldc < MAX(1, m))
		info = 13;
	නම් (ldb < MAX(1, nrowb))
		info = 10;
	නම් (lda < MAX(1, nrowa))
		info = 8;
	නම් (k < 0)
		info = 5;
	නම් (m < 0)
		info = 4;
	නම් (transb < 0)
		info = 3;
	නම් (transa < 0)
		info = 2;
	නම් (uplo < 0)
		info = 1;

	නම් (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ආපසු;
	}
#නැතහොත්

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	පූර්ණ_සංඛ්‍යා transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	නම් (order == CblasColMajor) {
		නම් (Uplo == CblasUpper) uplo = 0;
		නම් (Uplo == CblasLower) uplo = 1;

		නම් (TransA == CblasNoTrans)
			transa = 0;
		නම් (TransA == CblasTrans)
			transa = 1;

		නම් (TransA == CblasConjNoTrans)
			transa = 0;
		නම් (TransA == CblasConjTrans)
			transa = 1;

		නම් (TransB == CblasNoTrans)
			transb = 0;
		නම් (TransB == CblasTrans)
			transb = 1;

		නම් (TransB == CblasConjNoTrans)
			transb = 0;
		නම් (TransB == CblasConjTrans)
			transb = 1;

		a = (void *)A;
		b = (void *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		නම් (transa & 1) nrowa = k;
		nrowb = k;
		නම් (transb & 1)  nrowb = m;

		නම් (ldc < MAX(1, m))
			info = 13;
		නම් (ldb < MAX(1, nrowb))
			info = 10;
		නම් (lda < MAX(1, nrowa))
			info = 8;
		නම් (k < 0)
			info = 5;
		නම් (m < 0)
			info = 4;
		නම් (transb < 0)
			info = 3;
		නම් (transa < 0)
			info = 2;
		නම් (uplo < 0)
			info = 1;
	}

	නම් (order == CblasRowMajor) {

		a = (void *)B;
		b = (void *)A;

		lda = LDB;
		ldb = LDA;

		නම් (Uplo == CblasUpper) uplo = 0;
		නම් (Uplo == CblasLower) uplo = 1;

		නම් (TransB == CblasNoTrans)
			transa = 0;
		නම් (TransB == CblasTrans)
			transa = 1;

		නම් (TransB == CblasConjNoTrans)
			transa = 0;
		නම් (TransB == CblasConjTrans)
			transa = 1;

		නම් (TransA == CblasNoTrans)
			transb = 0;
		නම් (TransA == CblasTrans)
			transb = 1;

		නම් (TransA == CblasConjNoTrans)
			transb = 0;
		නම් (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		නම් (transa & 1) ncola = k;
		ncolb = k;

		නම් (transb & 1) {
			ncolb = m;
		}

		නම් (ldc < MAX(1,m))
			info = 13;
		නම් (ldb < MAX(1, ncolb))
			info = 8;
		නම් (lda < MAX(1, ncola))
			info = 10;
		නම් (k < 0)
			info = 5;
		නම් (m < 0)
			info = 4;
		නම් (transb < 0)
			info = 2;
		නම් (transa < 0)
			info = 3;
		නම් (uplo < 0)
			info = 1;
	}

	නම් (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, parimaa_na(ERROR_NAME));
		ආපසු;
	}

#endif
	පූර්ණ_සංඛ්‍යා buffer_size;
	blasint i, j;

#ifdef SMP
	පූර්ණ_සංඛ්‍යා nthreads;
#endif


#ifdef SMP
	static පූර්ණ_සංඛ්‍යා (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, පූර්ණ_සංඛ්‍යා) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	පූර්ණ_සංඛ්‍යා (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	නම් (m == 0)
		ආපසු;

	IDEBUG_START;

	const blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	නම් (uplo == 1) {
		සඳහා (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			නම් (transa & 1) {
				aa = a + lda * i;
			}
			නම් (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#නම් 0
			නම් (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			නම් (alpha == ZERO)
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

			නම් (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			නැතහොත්
				nthreads = num_cpu_avail(2);

			නම් (nthreads == 1) {
#endif

				නම් (!(transa & 1))
				(gemv[(පූර්ණ_සංඛ්‍යා)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				නැතහොත්
				(gemv[(පූර්ණ_සංඛ්‍යා)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} නැතහොත් {
				නම් (!(transa & 1))
				(gemv_thread[(පූර්ණ_සංඛ්‍යා)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				නැතහොත්
				(gemv_thread[(පූර්ණ_සංඛ්‍යා)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			STACK_FREE(buffer);
		}
	} නැතහොත් {

		සඳහා (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			නම් (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#නම් 0
			නම් (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, NULL, 0, NULL, 0);

			නම් (alpha == ZERO)
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

			නම් (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			නැතහොත්
				nthreads = num_cpu_avail(2);

			නම් (nthreads == 1) {
#endif

				නම් (!(transa & 1))
				(gemv[(පූර්ණ_සංඛ්‍යා)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				නැතහොත්
				(gemv[(පූර්ණ_සංඛ්‍යා)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} නැතහොත් {
				නම් (!(transa & 1))
				(gemv_thread[(පූර්ණ_සංඛ්‍යා)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				නැතහොත්
				(gemv_thread[(පූර්ණ_සංඛ්‍යා)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			STACK_FREE(buffer);
		}
	}

	IDEBUG_END;

	ආපසු;
}
