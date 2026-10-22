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

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن "common.h"

#اگر_تعریف_نشده COMPLEX
#تعریف_کن SMP_آستانه_MIN 65536.0
#اگر_تعریف_شده RNAME
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QGEMMTR"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DGEMMTR"
#وإلاإذا defined(BFLOAT16)
#تعریف_کن ERROR_NAME "SBGEMMTR"
#وإلا
#تعریف_کن ERROR_NAME "SGEMMTR"
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QGEMMT "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DGEMMT "
#وإلاإذا defined(BFLOAT16)
#تعریف_کن ERROR_NAME "SBGEMMT "
#وإلا
#تعریف_کن ERROR_NAME "SGEMMT "
#پایان_اگر
#پایان_اگر
#وإلا
#تعریف_کن SMP_آستانه_MIN 8192.0
#اگر_تعریف_شده RNAME
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XGEMMTR"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZGEMMTR"
#وإلا
#تعریف_کن ERROR_NAME "CGEMMTR"
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XGEMMT "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZGEMMT "
#وإلا
#تعریف_کن ERROR_NAME "CGEMMT "
#پایان_اگر
#پایان_اگر
#پایان_اگر

#اگر_تعریف_نشده GEMM_MULTITHREAD_آستانه
#تعریف_کن GEMM_MULTITHREAD_آستانه 4
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANSA, حرف *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	صحيح transa, transb, uplo;
	blasint info;

	حرف transA, transB, Uplo;
	blasint nrowa, nrowb;
#إذا !defined(COMPLEX)
	FLOAT alpha, beta;
#پایان_اگر

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#إذا defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#وإلا
	alpha = *Alpha;
	beta = *Beta;
#پایان_اگر

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

	إذا (transA == 'N')
		transa = 0;
	إذا (transA == 'T')
		transa = 1;
#اگر_تعریف_نشده COMPLEX
	إذا (transA == 'R')
		transa = 0;
	إذا (transA == 'C')
		transa = 1;
#وإلا
	إذا (transA == 'R')
		transa = 2;
	إذا (transA == 'C')
		transa = 3;
#پایان_اگر

	إذا (transB == 'N')
		transb = 0;
	إذا (transB == 'T')
		transb = 1;
#اگر_تعریف_نشده COMPLEX
	إذا (transB == 'R')
		transb = 0;
	إذا (transB == 'C')
		transb = 1;
#وإلا
	إذا (transB == 'R')
		transb = 2;
	إذا (transB == 'C')
		transb = 3;
#پایان_اگر

	إذا (Uplo == 'U')
		uplo = 0;
	إذا (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	إذا (transa & 1) nrowa = k;
	nrowb = k;
	إذا (transb & 1) nrowb = m;

	info = 0;

	إذا (ldc < MAX(1, m))
		info = 13;
	إذا (ldb < MAX(1, nrowb))
		info = 10;
	إذا (lda < MAX(1, nrowa))
		info = 8;
	إذا (k < 0)
		info = 5;
	إذا (m < 0)
		info = 4;
	إذا (transb < 0)
		info = 3;
	إذا (transa < 0)
		info = 2;
	إذا (uplo < 0)
		info = 1;

	إذا (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}
#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#اگر_تعریف_نشده COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint LDA,
	   فراغ *vb, blasint LDB, فراغ *vbeta, فراغ *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#پایان_اگر
	صحيح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	إذا (ترتیب == Cblas_ستون_Major) {
		إذا (Uplo == Cblas_فوقانی) uplo = 0;
		إذا (Uplo == Cblas_پایینی) uplo = 1;

		إذا (TransA == CblasNoTrans)
			transa = 0;
		إذا (TransA == CblasTrans)
			transa = 1;
#اگر_تعریف_نشده COMPLEX
		إذا (TransA == CblasConjNoTrans)
			transa = 0;
		إذا (TransA == CblasConjTrans)
			transa = 1;
#وإلا
		إذا (TransA == CblasConjNoTrans)
			transa = 2;
		إذا (TransA == CblasConjTrans)
			transa = 3;
#پایان_اگر
		إذا (TransB == CblasNoTrans)
			transb = 0;
		إذا (TransB == CblasTrans)
			transb = 1;
#اگر_تعریف_نشده COMPLEX
		إذا (TransB == CblasConjNoTrans)
			transb = 0;
		إذا (TransB == CblasConjTrans)
			transb = 1;
#وإلا
		إذا (TransB == CblasConjNoTrans)
			transb = 2;
		إذا (TransB == CblasConjTrans)
			transb = 3;
#پایان_اگر

		a = (فراغ *)A;
		b = (فراغ *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		إذا (transa & 1) nrowa = k;
		nrowb = k;
		إذا (transb & 1) nrowb = m;

		إذا (ldc < MAX(1, m))
			info = 13;
		إذا (ldb < MAX(1, nrowb))
			info = 10;
		إذا (lda < MAX(1, nrowa))
			info = 8;
		إذا (k < 0)
			info = 5;
		إذا (m < 0)
			info = 4;
		إذا (transb < 0)
			info = 3;
		إذا (transa < 0)
			info = 2;
		إذا (uplo < 0)
			info = 1;
	}

	إذا (ترتیب == Cblas_ردیف_Major) {

		a = (فراغ *)B;
		b = (فراغ *)A;

		lda = LDB;
		ldb = LDA;

		إذا (Uplo == Cblas_فوقانی) uplo = 1;
		إذا (Uplo == Cblas_پایینی) uplo = 0;

		إذا (TransB == CblasNoTrans)
			transa = 0;
		إذا (TransB == CblasTrans)
			transa = 1;
#اگر_تعریف_نشده COMPLEX
		إذا (TransB == CblasConjNoTrans)
			transa = 0;
		إذا (TransB == CblasConjTrans)
			transa = 1;
#وإلا
		إذا (TransB == CblasConjNoTrans)
			transa = 2;
		إذا (TransB == CblasConjTrans)
			transa = 3;
#پایان_اگر
		إذا (TransA == CblasNoTrans)
			transb = 0;
		إذا (TransA == CblasTrans)
			transb = 1;
#اگر_تعریف_نشده COMPLEX
		إذا (TransA == CblasConjNoTrans)
			transb = 0;
		إذا (TransA == CblasConjTrans)
			transb = 1;
#وإلا
		إذا (TransA == CblasConjNoTrans)
			transb = 2;
		إذا (TransA == CblasConjTrans)
			transb = 3;
#پایان_اگر

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		إذا (transa & 1) ncola = k;
		ncolb = k;
		إذا (transb & 1) ncolb = m;

		إذا (ldc < MAX(1,m))
			info = 13;
		إذا (ldb < MAX(1, ncolb))
			info = 8;
		إذا (lda < MAX(1, ncola))
			info = 10;
		إذا (k < 0)
			info = 5;
		إذا (m < 0)
			info = 4;
		إذا (transb < 0)
			info = 2;
		إذا (transa < 0)
			info = 3;
		إذا (uplo < 0)
			info = 1;
	}

	إذا (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}
#پایان_اگر
#اگر_تعریف_نشده COMPLEX
#تعریف_کن GEMMT_جدول_شاخص ((uplo << 2) | (transb << 1) | transa)
#وإلا
#تعریف_کن GEMMT_جدول_شاخص ((uplo << 4) | (transb << 2) | transa)
#پایان_اگر

	blas_arg_t آرگومانها;
	FLOAT *میانگیر, *sa, *sb;

#اگر_تعریف_شده SMP
	مزدوج MMK;
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
	صحيح نما = BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
	صحيح نما = BLAS_DOUBLE  | BLAS_REAL;
#وإلا
	صحيح نما = BLAS_SINGLE  | BLAS_REAL;
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
	صحيح نما = BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
	صحيح نما = BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
	صحيح نما = BLAS_SINGLE  | BLAS_COMPLEX;
#پایان_اگر
#پایان_اگر
#پایان_اگر

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ساكن صحيح (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_تعریف_نشده COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#وإلا
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#پایان_اگر
	};

	إذا (m == 0)
		أرجع;

	IDEBUG_شروع;

	FUNCTION_PROFILE_شروع();

	آرگومانها.m = m;
	آرگومانها.n = m;
	آرگومانها.k = k;

	آرگومانها.a = (فراغ *)a;
	آرگومانها.b = (فراغ *)b;
	آرگومانها.c = (فراغ *)c;

	آرگومانها.lda = lda;
	آرگومانها.ldb = ldb;
	آرگومانها.ldc = ldc;

#إذا defined(COMPLEX)
	آرگومانها.alpha = (فراغ *)alpha;
	آرگومانها.beta  = (فراغ *)beta;
#وإلا
	آرگومانها.alpha = (فراغ *)&alpha;
	آرگومانها.beta  = (فراغ *)&beta;
#پایان_اگر

#اگر_تعریف_شده SMP
	آرگومانها.common   = عدم;
	آرگومانها.nthreads = 1;
#پایان_اگر

	میانگیر = (FLOAT *)blas_حافظه_alloc(0);
	إذا (!میانگیر) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}

	sa = (FLOAT *)((BLASLONG)میانگیر + GEMM_جابه_جایی_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_جابه_جایی_B);

#اگر_تعریف_شده SMP
	MMK = (مزدوج)(m + 1) * (مزدوج)m * (مزدوج)k;
	إذا (MMK <= (SMP_آستانه_MIN * GEMM_MULTITHREAD_آستانه)) {
		آرگومانها.nthreads = 1;
	} وإلا {
		آرگومانها.nthreads = رقم_cpu_avail(3);
	}

	إذا (آرگومانها.nthreads == 1) {
#پایان_اگر

		(gemmt[GEMMT_جدول_شاخص])(&آرگومانها, عدم, عدم, sa, sb, 0);

#اگر_تعریف_شده SMP
	} وإلا {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		نما |= (uplo << BLAS_UPLO_لغزش);

		syrk_thread(نما, &آرگومانها, عدم, عدم, gemmt[GEMMT_جدول_شاخص], sa, sb, آرگومانها.nthreads);
	}
#پایان_اگر

	blas_حافظه_free(میانگیر);

	FUNCTION_PROFILE_انتها(COMPSIZE * COMPSIZE, آرگومانها.m * آرگومانها.k + آرگومانها.k * آرگومانها.n + آرگومانها.m * آرگومانها.n / 2,
			     آرگومانها.m * آرگومانها.n * آرگومانها.k);

	IDEBUG_انتها;

	أرجع;
}
