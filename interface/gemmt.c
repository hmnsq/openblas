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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ <stdlib.h>
#అన్తర్భావ "common.h"

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచించు SMP_పరిమితి_MIN 65536.0
#నిర్వచితమైతే RNAME
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGEMMTR"
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGEMMTR"
#లేదావేళ defined(BFLOAT16)
#నిర్వచించు ERROR_NAME "SBGEMMTR"
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGEMMTR"
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGEMMT "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGEMMT "
#లేదావేళ defined(BFLOAT16)
#నిర్వచించు ERROR_NAME "SBGEMMT "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGEMMT "
#అయితే_ముగింపు
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచించు SMP_పరిమితి_MIN 8192.0
#నిర్వచితమైతే RNAME
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGEMMTR"
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGEMMTR"
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGEMMTR"
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGEMMT "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGEMMT "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGEMMT "
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే GEMM_MULTITHREAD_పరిమితి
#నిర్వచించు GEMM_MULTITHREAD_పరిమితి 4
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, అక్షరం *TRANSA, అక్షరం *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	పూర్ణాంకం transa, transb, uplo;
	blasint info;

	అక్షరం transA, transB, Uplo;
	blasint nrowa, nrowb;
#ఒకవేళ !defined(COMPLEX)
	FLOAT alpha, beta;
#అయితే_ముగింపు

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#ఒకవేళ defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#లేకపోతే
	alpha = *Alpha;
	beta = *Beta;
#అయితే_ముగింపు

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

	ఒకవేళ (transA == 'N')
		transa = 0;
	ఒకవేళ (transA == 'T')
		transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
	ఒకవేళ (transA == 'R')
		transa = 0;
	ఒకవేళ (transA == 'C')
		transa = 1;
#లేకపోతే
	ఒకవేళ (transA == 'R')
		transa = 2;
	ఒకవేళ (transA == 'C')
		transa = 3;
#అయితే_ముగింపు

	ఒకవేళ (transB == 'N')
		transb = 0;
	ఒకవేళ (transB == 'T')
		transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
	ఒకవేళ (transB == 'R')
		transb = 0;
	ఒకవేళ (transB == 'C')
		transb = 1;
#లేకపోతే
	ఒకవేళ (transB == 'R')
		transb = 2;
	ఒకవేళ (transB == 'C')
		transb = 3;
#అయితే_ముగింపు

	ఒకవేళ (Uplo == 'U')
		uplo = 0;
	ఒకవేళ (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	ఒకవేళ (transa & 1) nrowa = k;
	nrowb = k;
	ఒకవేళ (transb & 1) nrowb = m;

	info = 0;

	ఒకవేళ (ldc < MAX(1, m))
		info = 13;
	ఒకవేళ (ldb < MAX(1, nrowb))
		info = 10;
	ఒకవేళ (lda < MAX(1, nrowa))
		info = 8;
	ఒకవేళ (k < 0)
		info = 5;
	ఒకవేళ (m < 0)
		info = 4;
	ఒకవేళ (transb < 0)
		info = 3;
	ఒకవేళ (transa < 0)
		info = 2;
	ఒకవేళ (uplo < 0)
		info = 1;

	ఒకవేళ (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
		తిరిగి;
	}
#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_UPLO Uplo,
	   గణన CBLAS_TRANSPOSE TransA, గణన CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#నిర్వచితం_కాకపోతే COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#లేకపోతే
	   శూన్య *valpha,
	   శూన్య *va, blasint LDA,
	   శూన్య *vb, blasint LDB, శూన్య *vbeta, శూన్య *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#అయితే_ముగింపు
	పూర్ణాంకం transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
		ఒకవేళ (Uplo == Cblas_పై) uplo = 0;
		ఒకవేళ (Uplo == Cblas_క్రింది) uplo = 1;

		ఒకవేళ (TransA == CblasNoTrans)
			transa = 0;
		ఒకవేళ (TransA == CblasTrans)
			transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
		ఒకవేళ (TransA == CblasConjNoTrans)
			transa = 0;
		ఒకవేళ (TransA == CblasConjTrans)
			transa = 1;
#లేకపోతే
		ఒకవేళ (TransA == CblasConjNoTrans)
			transa = 2;
		ఒకవేళ (TransA == CblasConjTrans)
			transa = 3;
#అయితే_ముగింపు
		ఒకవేళ (TransB == CblasNoTrans)
			transb = 0;
		ఒకవేళ (TransB == CblasTrans)
			transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
		ఒకవేళ (TransB == CblasConjNoTrans)
			transb = 0;
		ఒకవేళ (TransB == CblasConjTrans)
			transb = 1;
#లేకపోతే
		ఒకవేళ (TransB == CblasConjNoTrans)
			transb = 2;
		ఒకవేళ (TransB == CblasConjTrans)
			transb = 3;
#అయితే_ముగింపు

		a = (శూన్య *)A;
		b = (శూన్య *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		ఒకవేళ (transa & 1) nrowa = k;
		nrowb = k;
		ఒకవేళ (transb & 1) nrowb = m;

		ఒకవేళ (ldc < MAX(1, m))
			info = 13;
		ఒకవేళ (ldb < MAX(1, nrowb))
			info = 10;
		ఒకవేళ (lda < MAX(1, nrowa))
			info = 8;
		ఒకవేళ (k < 0)
			info = 5;
		ఒకవేళ (m < 0)
			info = 4;
		ఒకవేళ (transb < 0)
			info = 3;
		ఒకవేళ (transa < 0)
			info = 2;
		ఒకవేళ (uplo < 0)
			info = 1;
	}

	ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {

		a = (శూన్య *)B;
		b = (శూన్య *)A;

		lda = LDB;
		ldb = LDA;

		ఒకవేళ (Uplo == Cblas_పై) uplo = 1;
		ఒకవేళ (Uplo == Cblas_క్రింది) uplo = 0;

		ఒకవేళ (TransB == CblasNoTrans)
			transa = 0;
		ఒకవేళ (TransB == CblasTrans)
			transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
		ఒకవేళ (TransB == CblasConjNoTrans)
			transa = 0;
		ఒకవేళ (TransB == CblasConjTrans)
			transa = 1;
#లేకపోతే
		ఒకవేళ (TransB == CblasConjNoTrans)
			transa = 2;
		ఒకవేళ (TransB == CblasConjTrans)
			transa = 3;
#అయితే_ముగింపు
		ఒకవేళ (TransA == CblasNoTrans)
			transb = 0;
		ఒకవేళ (TransA == CblasTrans)
			transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
		ఒకవేళ (TransA == CblasConjNoTrans)
			transb = 0;
		ఒకవేళ (TransA == CblasConjTrans)
			transb = 1;
#లేకపోతే
		ఒకవేళ (TransA == CblasConjNoTrans)
			transb = 2;
		ఒకవేళ (TransA == CblasConjTrans)
			transb = 3;
#అయితే_ముగింపు

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		ఒకవేళ (transa & 1) ncola = k;
		ncolb = k;
		ఒకవేళ (transb & 1) ncolb = m;

		ఒకవేళ (ldc < MAX(1,m))
			info = 13;
		ఒకవేళ (ldb < MAX(1, ncolb))
			info = 8;
		ఒకవేళ (lda < MAX(1, ncola))
			info = 10;
		ఒకవేళ (k < 0)
			info = 5;
		ఒకవేళ (m < 0)
			info = 4;
		ఒకవేళ (transb < 0)
			info = 2;
		ఒకవేళ (transa < 0)
			info = 3;
		ఒకవేళ (uplo < 0)
			info = 1;
	}

	ఒకవేళ (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
		తిరిగి;
	}
#అయితే_ముగింపు
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచించు GEMMT_పట్టిక_సూచీ ((uplo << 2) | (transb << 1) | transa)
#లేకపోతే
#నిర్వచించు GEMMT_పట్టిక_సూచీ ((uplo << 4) | (transb << 2) | transa)
#అయితే_ముగింపు

	blas_arg_t వాదనలు;
	FLOAT *బఫరు, *sa, *sb;

#నిర్వచితమైతే SMP
	ద్విగుణ MMK;
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
	పూర్ణాంకం బహుళకం = BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
	పూర్ణాంకం బహుళకం = BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
	పూర్ణాంకం బహుళకం = BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
	పూర్ణాంకం బహుళకం = BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
	పూర్ణాంకం బహుళకం = BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
	పూర్ణాంకం బహుళకం = BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	స్థిర పూర్ణాంకం (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#నిర్వచితం_కాకపోతే COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#లేకపోతే
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#అయితే_ముగింపు
	};

	ఒకవేళ (m == 0)
		తిరిగి;

	IDEBUG_ప్రారంభం;

	FUNCTION_PROFILE_ప్రారంభం();

	వాదనలు.m = m;
	వాదనలు.n = m;
	వాదనలు.k = k;

	వాదనలు.a = (శూన్య *)a;
	వాదనలు.b = (శూన్య *)b;
	వాదనలు.c = (శూన్య *)c;

	వాదనలు.lda = lda;
	వాదనలు.ldb = ldb;
	వాదనలు.ldc = ldc;

#ఒకవేళ defined(COMPLEX)
	వాదనలు.alpha = (శూన్య *)alpha;
	వాదనలు.beta  = (శూన్య *)beta;
#లేకపోతే
	వాదనలు.alpha = (శూన్య *)&alpha;
	వాదనలు.beta  = (శూన్య *)&beta;
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
	వాదనలు.common   = ఖాళీ;
	వాదనలు.nthreads = 1;
#అయితే_ముగింపు

	బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(0);
	ఒకవేళ (!బఫరు) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
		తిరిగి;
	}

	sa = (FLOAT *)((BLASLONG)బఫరు + GEMM_స్థానభ్రంశం_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);

#నిర్వచితమైతే SMP
	MMK = (ద్విగుణ)(m + 1) * (ద్విగుణ)m * (ద్విగుణ)k;
	ఒకవేళ (MMK <= (SMP_పరిమితి_MIN * GEMM_MULTITHREAD_పరిమితి)) {
		వాదనలు.nthreads = 1;
	} లేకపోతే {
		వాదనలు.nthreads = అంకె_cpu_avail(3);
	}

	ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

		(gemmt[GEMMT_పట్టిక_సూచీ])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#నిర్వచితమైతే SMP
	} లేకపోతే {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		బహుళకం |= (uplo << BLAS_UPLO_జరుపు);

		syrk_thread(బహుళకం, &వాదనలు, ఖాళీ, ఖాళీ, gemmt[GEMMT_పట్టిక_సూచీ], sa, sb, వాదనలు.nthreads);
	}
#అయితే_ముగింపు

	blas_జ్ఞాపకం_free(బఫరు);

	FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE, వాదనలు.m * వాదనలు.k + వాదనలు.k * వాదనలు.n + వాదనలు.m * వాదనలు.n / 2,
			     వాదనలు.m * వాదనలు.n * వాదనలు.k);

	IDEBUG_చివర;

	తిరిగి;
}
