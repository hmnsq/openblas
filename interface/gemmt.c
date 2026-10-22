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

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#وصف_ڪريو SMP_حد_MIN 65536.0
#جيڪڏهن_وصف_ٿيل RNAME
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QGEMMTR"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DGEMMTR"
#وإلاإذا defined(BFLOAT16)
#وصف_ڪريو ERROR_NAME "SBGEMMTR"
#وإلا
#وصف_ڪريو ERROR_NAME "SGEMMTR"
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QGEMMT "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DGEMMT "
#وإلاإذا defined(BFLOAT16)
#وصف_ڪريو ERROR_NAME "SBGEMMT "
#وإلا
#وصف_ڪريو ERROR_NAME "SGEMMT "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#وإلا
#وصف_ڪريو SMP_حد_MIN 8192.0
#جيڪڏهن_وصف_ٿيل RNAME
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XGEMMTR"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZGEMMTR"
#وإلا
#وصف_ڪريو ERROR_NAME "CGEMMTR"
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XGEMMT "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZGEMMT "
#وإلا
#وصف_ڪريو ERROR_NAME "CGEMMT "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل GEMM_MULTITHREAD_حد
#وصف_ڪريو GEMM_MULTITHREAD_حد 4
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

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
#جيڪڏهن_پورو

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#إذا defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#وإلا
	alpha = *Alpha;
	beta = *Beta;
#جيڪڏهن_پورو

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
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	إذا (transA == 'R')
		transa = 0;
	إذا (transA == 'C')
		transa = 1;
#وإلا
	إذا (transA == 'R')
		transa = 2;
	إذا (transA == 'C')
		transa = 3;
#جيڪڏهن_پورو

	إذا (transB == 'N')
		transb = 0;
	إذا (transB == 'T')
		transb = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	إذا (transB == 'R')
		transb = 0;
	إذا (transB == 'C')
		transb = 1;
#وإلا
	إذا (transB == 'R')
		transb = 2;
	إذا (transB == 'C')
		transb = 3;
#جيڪڏهن_پورو

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

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
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
#جيڪڏهن_پورو
	صحيح transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	إذا (ترتيب == Cblas_کانو_Major) {
		إذا (Uplo == Cblas_مٿانهون) uplo = 0;
		إذا (Uplo == Cblas_هيٺاهون) uplo = 1;

		إذا (TransA == CblasNoTrans)
			transa = 0;
		إذا (TransA == CblasTrans)
			transa = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
		إذا (TransA == CblasConjNoTrans)
			transa = 0;
		إذا (TransA == CblasConjTrans)
			transa = 1;
#وإلا
		إذا (TransA == CblasConjNoTrans)
			transa = 2;
		إذا (TransA == CblasConjTrans)
			transa = 3;
#جيڪڏهن_پورو
		إذا (TransB == CblasNoTrans)
			transb = 0;
		إذا (TransB == CblasTrans)
			transb = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
		إذا (TransB == CblasConjNoTrans)
			transb = 0;
		إذا (TransB == CblasConjTrans)
			transb = 1;
#وإلا
		إذا (TransB == CblasConjNoTrans)
			transb = 2;
		إذا (TransB == CblasConjTrans)
			transb = 3;
#جيڪڏهن_پورو

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

	إذا (ترتيب == Cblas_صف_Major) {

		a = (فراغ *)B;
		b = (فراغ *)A;

		lda = LDB;
		ldb = LDA;

		إذا (Uplo == Cblas_مٿانهون) uplo = 1;
		إذا (Uplo == Cblas_هيٺاهون) uplo = 0;

		إذا (TransB == CblasNoTrans)
			transa = 0;
		إذا (TransB == CblasTrans)
			transa = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
		إذا (TransB == CblasConjNoTrans)
			transa = 0;
		إذا (TransB == CblasConjTrans)
			transa = 1;
#وإلا
		إذا (TransB == CblasConjNoTrans)
			transa = 2;
		إذا (TransB == CblasConjTrans)
			transa = 3;
#جيڪڏهن_پورو
		إذا (TransA == CblasNoTrans)
			transb = 0;
		إذا (TransA == CblasTrans)
			transb = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
		إذا (TransA == CblasConjNoTrans)
			transb = 0;
		إذا (TransA == CblasConjTrans)
			transb = 1;
#وإلا
		إذا (TransA == CblasConjNoTrans)
			transb = 2;
		إذا (TransA == CblasConjTrans)
			transb = 3;
#جيڪڏهن_پورو

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
#جيڪڏهن_پورو
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#وصف_ڪريو GEMMT_جدول_اشاريو ((uplo << 2) | (transb << 1) | transa)
#وإلا
#وصف_ڪريو GEMMT_جدول_اشاريو ((uplo << 4) | (transb << 2) | transa)
#جيڪڏهن_پورو

	blas_arg_t دليلون;
	FLOAT *بفر, *sa, *sb;

#جيڪڏهن_وصف_ٿيل SMP
	مزدوج MMK;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
	صحيح گهڻائي = BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
	صحيح گهڻائي = BLAS_DOUBLE  | BLAS_REAL;
#وإلا
	صحيح گهڻائي = BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
	صحيح گهڻائي = BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
	صحيح گهڻائي = BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
	صحيح گهڻائي = BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	ساكن صحيح (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
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
#جيڪڏهن_پورو
	};

	إذا (m == 0)
		أرجع;

	IDEBUG_شروعات;

	FUNCTION_PROFILE_شروعات();

	دليلون.m = m;
	دليلون.n = m;
	دليلون.k = k;

	دليلون.a = (فراغ *)a;
	دليلون.b = (فراغ *)b;
	دليلون.c = (فراغ *)c;

	دليلون.lda = lda;
	دليلون.ldb = ldb;
	دليلون.ldc = ldc;

#إذا defined(COMPLEX)
	دليلون.alpha = (فراغ *)alpha;
	دليلون.beta  = (فراغ *)beta;
#وإلا
	دليلون.alpha = (فراغ *)&alpha;
	دليلون.beta  = (فراغ *)&beta;
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
	دليلون.common   = عدم;
	دليلون.nthreads = 1;
#جيڪڏهن_پورو

	بفر = (FLOAT *)blas_يادگيري_alloc(0);
	إذا (!بفر) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, حجم(ERROR_NAME));
		أرجع;
	}

	sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);

#جيڪڏهن_وصف_ٿيل SMP
	MMK = (مزدوج)(m + 1) * (مزدوج)m * (مزدوج)k;
	إذا (MMK <= (SMP_حد_MIN * GEMM_MULTITHREAD_حد)) {
		دليلون.nthreads = 1;
	} وإلا {
		دليلون.nthreads = انگ_cpu_avail(3);
	}

	إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

		(gemmt[GEMMT_جدول_اشاريو])(&دليلون, عدم, عدم, sa, sb, 0);

#جيڪڏهن_وصف_ٿيل SMP
	} وإلا {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		گهڻائي |= (uplo << BLAS_UPLO_سرڪايو);

		syrk_thread(گهڻائي, &دليلون, عدم, عدم, gemmt[GEMMT_جدول_اشاريو], sa, sb, دليلون.nthreads);
	}
#جيڪڏهن_پورو

	blas_يادگيري_free(بفر);

	FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE, دليلون.m * دليلون.k + دليلون.k * دليلون.n + دليلون.m * دليلون.n / 2,
			     دليلون.m * دليلون.n * دليلون.k);

	IDEBUG_پڇاڙي;

	أرجع;
}
