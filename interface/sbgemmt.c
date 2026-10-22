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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <stdlib.h>
#അന്തര്ഭാവ "common.h"

#define SMP_THRESHOLD_MIN 65536.0
#define ERROR_NAME "SBGEMMT "

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, അക്ഷരം *TRANSA, അക്ഷരം *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	പൂർണ്ണസംഖ്യ transa, transb, uplo;
	blasint info;

	അക്ഷരം transA, transB, Uplo;
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

	എങ്കിൽ (transA == 'N')
		transa = 0;
	എങ്കിൽ (transA == 'T')
		transa = 1;

	എങ്കിൽ (transA == 'R')
		transa = 0;
	എങ്കിൽ (transA == 'C')
		transa = 1;

	എങ്കിൽ (transB == 'N')
		transb = 0;
	എങ്കിൽ (transB == 'T')
		transb = 1;

	എങ്കിൽ (transB == 'R')
		transb = 0;
	എങ്കിൽ (transB == 'C')
		transb = 1;

	എങ്കിൽ (Uplo == 'U')
		uplo = 0;
	എങ്കിൽ (Uplo == 'L')
		uplo = 1;
	nrowa = m;
	എങ്കിൽ (transa & 1) nrowa = k;
	nrowb = k;
	എങ്കിൽ (transb & 1) nrowb = m;

	info = 0;

	എങ്കിൽ (ldc < MAX(1, m))
		info = 13;
	എങ്കിൽ (ldb < MAX(1, nrowb))
		info = 10;
	എങ്കിൽ (lda < MAX(1, nrowa))
		info = 8;
	എങ്കിൽ (k < 0)
		info = 5;
	എങ്കിൽ (m < 0)
		info = 4;
	എങ്കിൽ (transb < 0)
		info = 3;
	എങ്കിൽ (transa < 0)
		info = 2;
	എങ്കിൽ (uplo < 0)
		info = 1;

	എങ്കിൽ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
		തിരികെ;
	}
#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA, എണ്ണൽ CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	പൂർണ്ണസംഖ്യ transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *buffer;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	എങ്കിൽ (order == CblasColMajor) {
		എങ്കിൽ (Uplo == CblasUpper) uplo = 0;
		എങ്കിൽ (Uplo == CblasLower) uplo = 1;

		എങ്കിൽ (TransA == CblasNoTrans)
			transa = 0;
		എങ്കിൽ (TransA == CblasTrans)
			transa = 1;

		എങ്കിൽ (TransA == CblasConjNoTrans)
			transa = 0;
		എങ്കിൽ (TransA == CblasConjTrans)
			transa = 1;

		എങ്കിൽ (TransB == CblasNoTrans)
			transb = 0;
		എങ്കിൽ (TransB == CblasTrans)
			transb = 1;

		എങ്കിൽ (TransB == CblasConjNoTrans)
			transb = 0;
		എങ്കിൽ (TransB == CblasConjTrans)
			transb = 1;

		a = (ശൂന്യം *)A;
		b = (ശൂന്യം *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		എങ്കിൽ (transa & 1) nrowa = k;
		nrowb = k;
		എങ്കിൽ (transb & 1)  nrowb = m;

		എങ്കിൽ (ldc < MAX(1, m))
			info = 13;
		എങ്കിൽ (ldb < MAX(1, nrowb))
			info = 10;
		എങ്കിൽ (lda < MAX(1, nrowa))
			info = 8;
		എങ്കിൽ (k < 0)
			info = 5;
		എങ്കിൽ (m < 0)
			info = 4;
		എങ്കിൽ (transb < 0)
			info = 3;
		എങ്കിൽ (transa < 0)
			info = 2;
		എങ്കിൽ (uplo < 0)
			info = 1;
	}

	എങ്കിൽ (order == CblasRowMajor) {

		a = (ശൂന്യം *)B;
		b = (ശൂന്യം *)A;

		lda = LDB;
		ldb = LDA;

		എങ്കിൽ (Uplo == CblasUpper) uplo = 0;
		എങ്കിൽ (Uplo == CblasLower) uplo = 1;

		എങ്കിൽ (TransB == CblasNoTrans)
			transa = 0;
		എങ്കിൽ (TransB == CblasTrans)
			transa = 1;

		എങ്കിൽ (TransB == CblasConjNoTrans)
			transa = 0;
		എങ്കിൽ (TransB == CblasConjTrans)
			transa = 1;

		എങ്കിൽ (TransA == CblasNoTrans)
			transb = 0;
		എങ്കിൽ (TransA == CblasTrans)
			transb = 1;

		എങ്കിൽ (TransA == CblasConjNoTrans)
			transb = 0;
		എങ്കിൽ (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		എങ്കിൽ (transa & 1) ncola = k;
		ncolb = k;

		എങ്കിൽ (transb & 1) {
			ncolb = m;
		}

		എങ്കിൽ (ldc < MAX(1,m))
			info = 13;
		എങ്കിൽ (ldb < MAX(1, ncolb))
			info = 8;
		എങ്കിൽ (lda < MAX(1, ncola))
			info = 10;
		എങ്കിൽ (k < 0)
			info = 5;
		എങ്കിൽ (m < 0)
			info = 4;
		എങ്കിൽ (transb < 0)
			info = 2;
		എങ്കിൽ (transa < 0)
			info = 3;
		എങ്കിൽ (uplo < 0)
			info = 1;
	}

	എങ്കിൽ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
		തിരികെ;
	}

#endif
	പൂർണ്ണസംഖ്യ buffer_size;
	blasint i, j;

#ifdef SMP
	പൂർണ്ണസംഖ്യ nthreads;
#endif


#ifdef SMP
	സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, പൂർണ്ണസംഖ്യ) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#endif
	പൂർണ്ണസംഖ്യ (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	എങ്കിൽ (m == 0)
		തിരികെ;

	IDEBUG_START;

	സ്ഥിരം blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	എങ്കിൽ (uplo == 1) {
		വേണ്ടി (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			എങ്കിൽ (transa & 1) {
				aa = a + lda * i;
			}
			എങ്കിൽ (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#എങ്കിൽ 0
			എങ്കിൽ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ശൂന്യ, 0, ശൂന്യ, 0);

			എങ്കിൽ (alpha == ZERO)
				തുടരുക;
#endif

			IDEBUG_START;

			buffer_size = j + k + 128 / വലുപ്പം(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / വലുപ്പം(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			അട്ടി_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			എങ്കിൽ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			അല്ലെങ്കിൽ
				nthreads = അക്കം_cpu_avail(2);

			എങ്കിൽ (nthreads == 1) {
#endif

				എങ്കിൽ (!(transa & 1))
				(gemv[(പൂർണ്ണസംഖ്യ)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				അല്ലെങ്കിൽ
				(gemv[(പൂർണ്ണസംഖ്യ)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#ifdef SMP
			} അല്ലെങ്കിൽ {
				എങ്കിൽ (!(transa & 1))
				(gemv_thread[(പൂർണ്ണസംഖ്യ)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				അല്ലെങ്കിൽ
				(gemv_thread[(പൂർണ്ണസംഖ്യ)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#endif

			അട്ടി_FREE(buffer);
		}
	} അല്ലെങ്കിൽ {

		വേണ്ടി (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			എങ്കിൽ (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#എങ്കിൽ 0
			എങ്കിൽ (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, ശൂന്യ, 0, ശൂന്യ, 0);

			എങ്കിൽ (alpha == ZERO)
				തുടരുക;
#endif
			IDEBUG_START;

			buffer_size = j + k + 128 / വലുപ്പം(FLOAT);
#ifdef WINDOWS_ABI
			buffer_size += 160 / വലുപ്പം(FLOAT);
#endif
			// for alignment
			buffer_size = (buffer_size + 3) & ~3;
			അട്ടി_ALLOC(buffer_size, IFLOAT, buffer);

#ifdef SMP

			എങ്കിൽ (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			അല്ലെങ്കിൽ
				nthreads = അക്കം_cpu_avail(2);

			എങ്കിൽ (nthreads == 1) {
#endif

				എങ്കിൽ (!(transa & 1))
				(gemv[(പൂർണ്ണസംഖ്യ)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				അല്ലെങ്കിൽ
				(gemv[(പൂർണ്ണസംഖ്യ)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#ifdef SMP
			} അല്ലെങ്കിൽ {
				എങ്കിൽ (!(transa & 1))
				(gemv_thread[(പൂർണ്ണസംഖ്യ)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				അല്ലെങ്കിൽ
				(gemv_thread[(പൂർണ്ണസംഖ്യ)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#endif

			അട്ടി_FREE(buffer);
		}
	}

	IDEBUG_END;

	തിരികെ;
}
