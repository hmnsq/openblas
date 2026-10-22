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

#که_نه_تعریف_شوی COMPLEX
#تعریف_کړه SMP_پوله_MIN 65536.0
#که_تعریف_شوی RNAME
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QGEMMTR"
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DGEMMTR"
#وإلاإذا defined(BFLOAT16)
#تعریف_کړه ERROR_NAME "SBGEMMTR"
#وإلا
#تعریف_کړه ERROR_NAME "SGEMMTR"
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QGEMMT "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DGEMMT "
#وإلاإذا defined(BFLOAT16)
#تعریف_کړه ERROR_NAME "SBGEMMT "
#وإلا
#تعریف_کړه ERROR_NAME "SGEMMT "
#د_که_پای
#د_که_پای
#وإلا
#تعریف_کړه SMP_پوله_MIN 8192.0
#که_تعریف_شوی RNAME
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XGEMMTR"
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZGEMMTR"
#وإلا
#تعریف_کړه ERROR_NAME "CGEMMTR"
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XGEMMT "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZGEMMT "
#وإلا
#تعریف_کړه ERROR_NAME "CGEMMT "
#د_که_پای
#د_که_پای
#د_که_پای

#که_نه_تعریف_شوی GEMM_MULTITHREAD_پوله
#تعریف_کړه GEMM_MULTITHREAD_پوله 4
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

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
#د_که_پای

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#إذا defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#وإلا
	alpha = *Alpha;
	beta = *Beta;
#د_که_پای

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
#که_نه_تعریف_شوی COMPLEX
	إذا (transA == 'R')
		transa = 0;
	إذا (transA == 'C')
		transa = 1;
#وإلا
	إذا (transA == 'R')
		transa = 2;
	إذا (transA == 'C')
		transa = 3;
#د_که_پای

	إذا (transB == 'N')
		transb = 0;
	إذا (transB == 'T')
		transb = 1;
#که_نه_تعریف_شوی COMPLEX
	إذا (transB == 'R')
		transb = 0;
	إذا (transB == 'C')
		transb = 1;
#وإلا
	إذا (transB == 'R')
		transb = 2;
	إذا (transB == 'C')
		transb = 3;
#د_که_پای

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
#که_نه_تعریف_شوی COMPLEX
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
#د_که_پای
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
		إذا (Uplo == Cblas_پاسنی) uplo = 0;
		إذا (Uplo == Cblas_ښکتنی) uplo = 1;

		إذا (TransA == CblasNoTrans)
			transa = 0;
		إذا (TransA == CblasTrans)
			transa = 1;
#که_نه_تعریف_شوی COMPLEX
		إذا (TransA == CblasConjNoTrans)
			transa = 0;
		إذا (TransA == CblasConjTrans)
			transa = 1;
#وإلا
		إذا (TransA == CblasConjNoTrans)
			transa = 2;
		إذا (TransA == CblasConjTrans)
			transa = 3;
#د_که_پای
		إذا (TransB == CblasNoTrans)
			transb = 0;
		إذا (TransB == CblasTrans)
			transb = 1;
#که_نه_تعریف_شوی COMPLEX
		إذا (TransB == CblasConjNoTrans)
			transb = 0;
		إذا (TransB == CblasConjTrans)
			transb = 1;
#وإلا
		إذا (TransB == CblasConjNoTrans)
			transb = 2;
		إذا (TransB == CblasConjTrans)
			transb = 3;
#د_که_پای

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

	إذا (ترتیب == Cblas_لیکه_Major) {

		a = (فراغ *)B;
		b = (فراغ *)A;

		lda = LDB;
		ldb = LDA;

		إذا (Uplo == Cblas_پاسنی) uplo = 1;
		إذا (Uplo == Cblas_ښکتنی) uplo = 0;

		إذا (TransB == CblasNoTrans)
			transa = 0;
		إذا (TransB == CblasTrans)
			transa = 1;
#که_نه_تعریف_شوی COMPLEX
		إذا (TransB == CblasConjNoTrans)
			transa = 0;
		إذا (TransB == CblasConjTrans)
			transa = 1;
#وإلا
		إذا (TransB == CblasConjNoTrans)
			transa = 2;
		إذا (TransB == CblasConjTrans)
			transa = 3;
#د_که_پای
		إذا (TransA == CblasNoTrans)
			transb = 0;
		إذا (TransA == CblasTrans)
			transb = 1;
#که_نه_تعریف_شوی COMPLEX
		إذا (TransA == CblasConjNoTrans)
			transb = 0;
		إذا (TransA == CblasConjTrans)
			transb = 1;
#وإلا
		إذا (TransA == CblasConjNoTrans)
			transb = 2;
		إذا (TransA == CblasConjTrans)
			transb = 3;
#د_که_پای

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
#د_که_پای
#که_نه_تعریف_شوی COMPLEX
#تعریف_کړه GEMMT_جدول_شاخص ((uplo << 2) | (transb << 1) | transa)
#وإلا
#تعریف_کړه GEMMT_جدول_شاخص ((uplo << 4) | (transb << 2) | transa)
#د_که_پای

	blas_arg_t دلیلونه;
	FLOAT *بفر, *sa, *sb;

#که_تعریف_شوی SMP
	مزدوج MMK;
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
	صحيح ډېرپېښ = BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
	صحيح ډېرپېښ = BLAS_DOUBLE  | BLAS_REAL;
#وإلا
	صحيح ډېرپېښ = BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
	صحيح ډېرپېښ = BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
	صحيح ډېرپېښ = BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
	صحيح ډېرپېښ = BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ساكن صحيح (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#که_نه_تعریف_شوی COMPLEX
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
#د_که_پای
	};

	إذا (m == 0)
		أرجع;

	IDEBUG_پیل;

	FUNCTION_PROFILE_پیل();

	دلیلونه.m = m;
	دلیلونه.n = m;
	دلیلونه.k = k;

	دلیلونه.a = (فراغ *)a;
	دلیلونه.b = (فراغ *)b;
	دلیلونه.c = (فراغ *)c;

	دلیلونه.lda = lda;
	دلیلونه.ldb = ldb;
	دلیلونه.ldc = ldc;

#إذا defined(COMPLEX)
	دلیلونه.alpha = (فراغ *)alpha;
	دلیلونه.beta  = (فراغ *)beta;
#وإلا
	دلیلونه.alpha = (فراغ *)&alpha;
	دلیلونه.beta  = (فراغ *)&beta;
#د_که_پای

#که_تعریف_شوی SMP
	دلیلونه.common   = عدم;
	دلیلونه.nthreads = 1;
#د_که_پای

	بفر = (FLOAT *)blas_حافظه_alloc(0);
	إذا (!بفر) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}

	sa = (FLOAT *)((BLASLONG)بفر + GEMM_دننه_واټن_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_دننه_واټن_B);

#که_تعریف_شوی SMP
	MMK = (مزدوج)(m + 1) * (مزدوج)m * (مزدوج)k;
	إذا (MMK <= (SMP_پوله_MIN * GEMM_MULTITHREAD_پوله)) {
		دلیلونه.nthreads = 1;
	} وإلا {
		دلیلونه.nthreads = عدد_cpu_avail(3);
	}

	إذا (دلیلونه.nthreads == 1) {
#د_که_پای

		(gemmt[GEMMT_جدول_شاخص])(&دلیلونه, عدم, عدم, sa, sb, 0);

#که_تعریف_شوی SMP
	} وإلا {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		ډېرپېښ |= (uplo << BLAS_UPLO_ښويدنه);

		syrk_thread(ډېرپېښ, &دلیلونه, عدم, عدم, gemmt[GEMMT_جدول_شاخص], sa, sb, دلیلونه.nthreads);
	}
#د_که_پای

	blas_حافظه_free(بفر);

	FUNCTION_PROFILE_پای_ټکی(COMPSIZE * COMPSIZE, دلیلونه.m * دلیلونه.k + دلیلونه.k * دلیلونه.n + دلیلونه.m * دلیلونه.n / 2,
			     دلیلونه.m * دلیلونه.n * دلیلونه.k);

	IDEBUG_پای_ټکی;

	أرجع;
}
