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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <stdlib.h>
#અન્તર્ભાવ "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, અક્ષર *TRANSA, અક્ષર *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	પૂર્ણાંક transa, transb, uplo;
	blasint info;

	અક્ષર transA, transB, Uplo;
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

	જો (transA == 'N')
		transa = 0;
	જો (transA == 'T')
		transa = 1;

	જો (transA == 'R')
		transa = 0;
	જો (transA == 'C')
		transa = 1;

	જો (transB == 'N')
		transb = 0;
	જો (transB == 'T')
		transb = 1;

	જો (transB == 'R')
		transb = 0;
	જો (transB == 'C')
		transb = 1;

	જો (Uplo == 'U')
		uplo = 0;
	જો (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	જો (transa & 1) nrowa = k;
	nrowb = k;
	જો (transb & 1) nrowb = m;

	info = 0;

	જો (ldc < MAX(1, m))
		info = 13;
	જો (ldb < MAX(1, nrowb))
		info = 10;
	જો (lda < MAX(1, nrowa))
		info = 8;
	જો (k < 0)
		info = 5;
	જો (m < 0)
		info = 4;
	જો (transb < 0)
		info = 3;
	જો (transa < 0)
		info = 2;
	જો (uplo < 0)
		info = 1;

	જો (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, કદ(ERROR_NAME));
		પાછા;
	}
#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_UPLO Uplo,
	   ગણના CBLAS_TRANSPOSE TransA, ગણના CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	પૂર્ણાંક transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	જો (order == CblasColMajor) {
		જો (Uplo == CblasUpper) uplo = 0;
		જો (Uplo == CblasLower) uplo = 1;

		જો (TransA == CblasNoTrans)
			transa = 0;
		જો (TransA == CblasTrans)
			transa = 1;

		જો (TransA == CblasConjNoTrans)
			transa = 0;
		જો (TransA == CblasConjTrans)
			transa = 1;

		જો (TransB == CblasNoTrans)
			transb = 0;
		જો (TransB == CblasTrans)
			transb = 1;

		જો (TransB == CblasConjNoTrans)
			transb = 0;
		જો (TransB == CblasConjTrans)
			transb = 1;

		a = (શૂન્ય *)A;
		b = (શૂન્ય *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		જો (transa & 1) nrowa = k;
		nrowb = k;
		જો (transb & 1)  nrowb = m;

		જો (ldc < MAX(1, m))
			info = 13;
		જો (ldb < MAX(1, nrowb))
			info = 10;
		જો (lda < MAX(1, nrowa))
			info = 8;
		જો (k < 0)
			info = 5;
		જો (m < 0)
			info = 4;
		જો (transb < 0)
			info = 3;
		જો (transa < 0)
			info = 2;
		જો (uplo < 0)
			info = 1;
	}

	જો (order == CblasRowMajor) {

		a = (શૂન્ય *)B;
		b = (શૂન્ય *)A;

		lda = LDB;
		ldb = LDA;

		જો (Uplo == CblasUpper) uplo = 0;
		જો (Uplo == CblasLower) uplo = 1;

		જો (TransB == CblasNoTrans)
			transa = 0;
		જો (TransB == CblasTrans)
			transa = 1;

		જો (TransB == CblasConjNoTrans)
			transa = 0;
		જો (TransB == CblasConjTrans)
			transa = 1;

		જો (TransA == CblasNoTrans)
			transb = 0;
		જો (TransA == CblasTrans)
			transb = 1;

		જો (TransA == CblasConjNoTrans)
			transb = 0;
		જો (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		જો (transa & 1) ncola = k;
		ncolb = k;

		જો (transb & 1) {
			ncolb = m;
		}

		જો (ldc < MAX(1,m))
			info = 13;
		જો (ldb < MAX(1, ncolb))
			info = 8;
		જો (lda < MAX(1, ncola))
			info = 10;
		જો (k < 0)
			info = 5;
		જો (m < 0)
			info = 4;
		જો (transb < 0)
			info = 2;
		જો (transa < 0)
			info = 3;
		જો (uplo < 0)
			info = 1;
	}

	જો (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, કદ(ERROR_NAME));
		પાછા;
	}

#endif
	પૂર્ણાંક buffer_size;
	blasint i, j;

#ifdef SMP
	પૂર્ણાંક nthreads;
#endif


#ifdef SMP
	સ્થિર પૂર્ણાંક (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, પૂર્ણાંક) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	પૂર્ણાંક (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	જો (m == 0)
		પાછા;

	IDEBUG_START;

	અચલ blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	જો (uplo == 1) {
		માટે (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			જો (transa & 1) {
				aa = a + lda * i;
			}
			જો (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#જો 0
			જો (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ખાલી, 0, ખાલી, 0);

			જો (alpha == ZERO)
				ચાલુ;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / કદ(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / કદ(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			ઢગલો_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			જો (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			નહીંતર
				nthreads = અંક_cpu_avail(2);

			જો (nthreads == 1) {
#endif

				જો (!(transa & 1))
				(gemv[(પૂર્ણાંક)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				નહીંતર
				(gemv[(પૂર્ણાંક)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} નહીંતર {
				જો (!(transa & 1))
				(gemv_thread[(પૂર્ણાંક)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				નહીંતર
				(gemv_thread[(પૂર્ણાંક)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			ઢગલો_FREE(buffer);
		}
	} નહીંતર {

		માટે (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			જો (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#જો 0
			જો (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ખાલી, 0, ખાલી, 0);

			જો (alpha == ZERO)
				ચાલુ;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / કદ(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / કદ(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			ઢગલો_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			જો (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			નહીંતર
				nthreads = અંક_cpu_avail(2);

			જો (nthreads == 1) {
#endif

				જો (!(transa & 1))
				(gemv[(પૂર્ણાંક)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				નહીંતર
				(gemv[(પૂર્ણાંક)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} નહીંતર {
				જો (!(transa & 1))
				(gemv_thread[(પૂર્ણાંક)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				નહીંતર
				(gemv_thread[(પૂર્ણાંક)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			ઢગલો_FREE(buffer);
		}
	}

	IDEBUG_END;

	પાછા;
}
