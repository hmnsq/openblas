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

#تعریف_کرو SMP_حد_MIN 65536.0
#تعریف_کرو ERROR_NAME "SBGEMMT "

#اگر_غیر_تعریف_شدہ GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *UPLO, حرف *TRANSA, حرف *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	صحیح transa, transb, uplo;
	blasint info;

	حرف transA, transB, Uplo;
	blasint nrowa, nrowb;
	IFLOAT *بفر;
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
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		واپس;
	}
#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_UPLO Uplo,
	   شمار CBLAS_TRANSPOSE TransA, شمار CBLAS_TRANSPOSE TransB, blasint m,
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
	XFLOAT *بفر;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	اگر (ترتیب == Cblas_ستون_Major) {
		اگر (Uplo == Cblas_بالائی) uplo = 0;
		اگر (Uplo == Cblas_زیریں) uplo = 1;

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

		a = (خالی *)A;
		b = (خالی *)B;
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

	اگر (ترتیب == Cblas_صف_Major) {

		a = (خالی *)B;
		b = (خالی *)A;

		lda = LDB;
		ldb = LDA;

		اگر (Uplo == Cblas_بالائی) uplo = 0;
		اگر (Uplo == Cblas_زیریں) uplo = 1;

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
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		واپس;
	}

#اگر_ختم
	صحیح بفر_سائز;
	blasint i, j;

#اگر_تعریف_شدہ SMP
	صحیح nthreads;
#اگر_ختم


#اگر_تعریف_شدہ SMP
	ساکن صحیح (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, صحیح) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#اگر_ختم
	صحیح (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	اگر (m == 0)
		واپس;

	IDEBUG_آغاز;

	ثابت blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	اگر (uplo == 1) {
		برائے (i = 0; i < m; i++) {
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
				SCAL_K(l, 0, 0, beta, cc, 1, تہی, 0, تہی, 0);

			اگر (alpha == ZERO)
				جاری;
#اگر_ختم

			IDEBUG_آغاز;

			بفر_سائز = j + k + 128 / حجم(FLOAT);
#اگر_تعریف_شدہ WINDOWS_ABI
			بفر_سائز += 160 / حجم(FLOAT);
#اگر_ختم
			// for alignment
			بفر_سائز = (بفر_سائز + 3) & ~3;
			ڈھیر_ALLOC(بفر_سائز, IFLOAT, بفر);

#اگر_تعریف_شدہ SMP

			اگر (1L * j * k < 2304L * GEMM_MULTITHREAD_حد)
				nthreads = 1;
			ورنہ
				nthreads = عدد_cpu_avail(2);

			اگر (nthreads == 1) {
#اگر_ختم

				اگر (!(transa & 1))
				(gemv[(صحیح)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				ورنہ
				(gemv[(صحیح)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#اگر_تعریف_شدہ SMP
			} ورنہ {
				اگر (!(transa & 1))
				(gemv_thread[(صحیح)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				ورنہ
				(gemv_thread[(صحیح)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#اگر_ختم

			ڈھیر_FREE(بفر);
		}
	} ورنہ {

		برائے (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			اگر (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#اگر 0
			اگر (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, تہی, 0, تہی, 0);

			اگر (alpha == ZERO)
				جاری;
#اگر_ختم
			IDEBUG_آغاز;

			بفر_سائز = j + k + 128 / حجم(FLOAT);
#اگر_تعریف_شدہ WINDOWS_ABI
			بفر_سائز += 160 / حجم(FLOAT);
#اگر_ختم
			// for alignment
			بفر_سائز = (بفر_سائز + 3) & ~3;
			ڈھیر_ALLOC(بفر_سائز, IFLOAT, بفر);

#اگر_تعریف_شدہ SMP

			اگر (1L * j * k < 2304L * GEMM_MULTITHREAD_حد)
				nthreads = 1;
			ورنہ
				nthreads = عدد_cpu_avail(2);

			اگر (nthreads == 1) {
#اگر_ختم

				اگر (!(transa & 1))
				(gemv[(صحیح)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				ورنہ
				(gemv[(صحیح)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#اگر_تعریف_شدہ SMP
			} ورنہ {
				اگر (!(transa & 1))
				(gemv_thread[(صحیح)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				ورنہ
				(gemv_thread[(صحیح)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#اگر_ختم

			ڈھیر_FREE(بفر);
		}
	}

	IDEBUG_انتہا;

	واپس;
}
