/*********************************************************************/
/* Copyright 2022, The OpenBLAS Project.                             */
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

#اگر_غیر_تعریف_شدہ COMPLEX
#تعریف_کرو SMP_حد_MIN 65536.0
#اگر_تعریف_شدہ RNAME
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMMTR"
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMMTR"
#وگرنہ_اگر defined(BFLOAT16)
#تعریف_کرو ERROR_NAME "SBGEMMTR"
#ورنہ
#تعریف_کرو ERROR_NAME "SGEMMTR"
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMMT "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMMT "
#وگرنہ_اگر defined(BFLOAT16)
#تعریف_کرو ERROR_NAME "SBGEMMT "
#ورنہ
#تعریف_کرو ERROR_NAME "SGEMMT "
#اگر_ختم
#اگر_ختم
#ورنہ
#تعریف_کرو SMP_حد_MIN 8192.0
#اگر_تعریف_شدہ RNAME
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMMTR"
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMMTR"
#ورنہ
#تعریف_کرو ERROR_NAME "CGEMMTR"
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMMT "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMMT "
#ورنہ
#تعریف_کرو ERROR_NAME "CGEMMT "
#اگر_ختم
#اگر_ختم
#اگر_ختم

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
#اگر !defined(COMPLEX)
	FLOAT alpha, beta;
#اگر_ختم

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#اگر defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#ورنہ
	alpha = *Alpha;
	beta = *Beta;
#اگر_ختم

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
#اگر_غیر_تعریف_شدہ COMPLEX
	اگر (transA == 'R')
		transa = 0;
	اگر (transA == 'C')
		transa = 1;
#ورنہ
	اگر (transA == 'R')
		transa = 2;
	اگر (transA == 'C')
		transa = 3;
#اگر_ختم

	اگر (transB == 'N')
		transb = 0;
	اگر (transB == 'T')
		transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
	اگر (transB == 'R')
		transb = 0;
	اگر (transB == 'C')
		transb = 1;
#ورنہ
	اگر (transB == 'R')
		transb = 2;
	اگر (transB == 'C')
		transb = 3;
#اگر_ختم

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
#اگر_غیر_تعریف_شدہ COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#ورنہ
	   خالی *valpha,
	   خالی *va, blasint LDA,
	   خالی *vb, blasint LDB, خالی *vbeta, خالی *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#اگر_ختم
	صحیح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

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
#اگر_غیر_تعریف_شدہ COMPLEX
		اگر (TransA == CblasConjNoTrans)
			transa = 0;
		اگر (TransA == CblasConjTrans)
			transa = 1;
#ورنہ
		اگر (TransA == CblasConjNoTrans)
			transa = 2;
		اگر (TransA == CblasConjTrans)
			transa = 3;
#اگر_ختم
		اگر (TransB == CblasNoTrans)
			transb = 0;
		اگر (TransB == CblasTrans)
			transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
		اگر (TransB == CblasConjNoTrans)
			transb = 0;
		اگر (TransB == CblasConjTrans)
			transb = 1;
#ورنہ
		اگر (TransB == CblasConjNoTrans)
			transb = 2;
		اگر (TransB == CblasConjTrans)
			transb = 3;
#اگر_ختم

		a = (خالی *)A;
		b = (خالی *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		اگر (transa & 1) nrowa = k;
		nrowb = k;
		اگر (transb & 1) nrowb = m;

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

		اگر (Uplo == Cblas_بالائی) uplo = 1;
		اگر (Uplo == Cblas_زیریں) uplo = 0;

		اگر (TransB == CblasNoTrans)
			transa = 0;
		اگر (TransB == CblasTrans)
			transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
		اگر (TransB == CblasConjNoTrans)
			transa = 0;
		اگر (TransB == CblasConjTrans)
			transa = 1;
#ورنہ
		اگر (TransB == CblasConjNoTrans)
			transa = 2;
		اگر (TransB == CblasConjTrans)
			transa = 3;
#اگر_ختم
		اگر (TransA == CblasNoTrans)
			transb = 0;
		اگر (TransA == CblasTrans)
			transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
		اگر (TransA == CblasConjNoTrans)
			transb = 0;
		اگر (TransA == CblasConjTrans)
			transb = 1;
#ورنہ
		اگر (TransA == CblasConjNoTrans)
			transb = 2;
		اگر (TransA == CblasConjTrans)
			transb = 3;
#اگر_ختم

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		اگر (transa & 1) ncola = k;
		ncolb = k;
		اگر (transb & 1) ncolb = m;

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
#اگر_غیر_تعریف_شدہ COMPLEX
#تعریف_کرو GEMMT_جدول_اشاریہ ((uplo << 2) | (transb << 1) | transa)
#ورنہ
#تعریف_کرو GEMMT_جدول_اشاریہ ((uplo << 4) | (transb << 2) | transa)
#اگر_ختم

	blas_arg_t دلائل;
	FLOAT *بفر, *sa, *sb;

#اگر_تعریف_شدہ SMP
	دوہرا MMK;
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
	صحیح منوال = BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
	صحیح منوال = BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
	صحیح منوال = BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
	صحیح منوال = BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
	صحیح منوال = BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
	صحیح منوال = BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ساکن صحیح (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_غیر_تعریف_شدہ COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#ورنہ
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#اگر_ختم
	};

	اگر (m == 0)
		واپس;

	IDEBUG_آغاز;

	FUNCTION_PROFILE_آغاز();

	دلائل.m = m;
	دلائل.n = m;
	دلائل.k = k;

	دلائل.a = (خالی *)a;
	دلائل.b = (خالی *)b;
	دلائل.c = (خالی *)c;

	دلائل.lda = lda;
	دلائل.ldb = ldb;
	دلائل.ldc = ldc;

#اگر defined(COMPLEX)
	دلائل.alpha = (خالی *)alpha;
	دلائل.beta  = (خالی *)beta;
#ورنہ
	دلائل.alpha = (خالی *)&alpha;
	دلائل.beta  = (خالی *)&beta;
#اگر_ختم

#اگر_تعریف_شدہ SMP
	دلائل.common   = تہی;
	دلائل.nthreads = 1;
#اگر_ختم

	بفر = (FLOAT *)blas_حافظہ_alloc(0);
	اگر (!بفر) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		واپس;
	}

	sa = (FLOAT *)((BLASLONG)بفر + GEMM_انحراف_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);

#اگر_تعریف_شدہ SMP
	MMK = (دوہرا)(m + 1) * (دوہرا)m * (دوہرا)k;
	اگر (MMK <= (SMP_حد_MIN * GEMM_MULTITHREAD_حد)) {
		دلائل.nthreads = 1;
	} ورنہ {
		دلائل.nthreads = عدد_cpu_avail(3);
	}

	اگر (دلائل.nthreads == 1) {
#اگر_ختم

		(gemmt[GEMMT_جدول_اشاریہ])(&دلائل, تہی, تہی, sa, sb, 0);

#اگر_تعریف_شدہ SMP
	} ورنہ {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		منوال |= (uplo << BLAS_UPLO_سرکاؤ);

		syrk_thread(منوال, &دلائل, تہی, تہی, gemmt[GEMMT_جدول_اشاریہ], sa, sb, دلائل.nthreads);
	}
#اگر_ختم

	blas_حافظہ_free(بفر);

	FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE, دلائل.m * دلائل.k + دلائل.k * دلائل.n + دلائل.m * دلائل.n / 2,
			     دلائل.m * دلائل.n * دلائل.k);

	IDEBUG_انتہا;

	واپس;
}
