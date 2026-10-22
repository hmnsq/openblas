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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव "common.h"

#यदि_अपरिभाषितम् COMPLEX
#परिभाषय SMP_देहली_MIN 65536.0
#यदि_परिभाषितम् RNAME
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGEMMTR"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGEMMTR"
#नोचेद्यदि defined(BFLOAT16)
#परिभाषय ERROR_NAME "SBGEMMTR"
#अन्यथा
#परिभाषय ERROR_NAME "SGEMMTR"
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGEMMT "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGEMMT "
#नोचेद्यदि defined(BFLOAT16)
#परिभाषय ERROR_NAME "SBGEMMT "
#अन्यथा
#परिभाषय ERROR_NAME "SGEMMT "
#यदि_अन्तः
#यदि_अन्तः
#अन्यथा
#परिभाषय SMP_देहली_MIN 8192.0
#यदि_परिभाषितम् RNAME
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMMTR"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMMTR"
#अन्यथा
#परिभाषय ERROR_NAME "CGEMMTR"
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMMT "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMMT "
#अन्यथा
#परिभाषय ERROR_NAME "CGEMMT "
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि_अपरिभाषितम् GEMM_MULTITHREAD_देहली
#परिभाषय GEMM_MULTITHREAD_देहली 4
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

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
#यदि !defined(COMPLEX)
	FLOAT alpha, beta;
#यदि_अन्तः

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#यदि defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#अन्यथा
	alpha = *Alpha;
	beta = *Beta;
#यदि_अन्तः

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
#यदि_अपरिभाषितम् COMPLEX
	यदि (transA == 'R')
		transa = 0;
	यदि (transA == 'C')
		transa = 1;
#अन्यथा
	यदि (transA == 'R')
		transa = 2;
	यदि (transA == 'C')
		transa = 3;
#यदि_अन्तः

	यदि (transB == 'N')
		transb = 0;
	यदि (transB == 'T')
		transb = 1;
#यदि_अपरिभाषितम् COMPLEX
	यदि (transB == 'R')
		transb = 0;
	यदि (transB == 'C')
		transb = 1;
#अन्यथा
	यदि (transB == 'R')
		transb = 2;
	यदि (transB == 'C')
		transb = 3;
#यदि_अन्तः

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

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint LDA,
	   शून्य *vb, blasint LDB, शून्य *vbeta, शून्य *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#यदि_अन्तः
	पूर्णाङ्क transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	यदि (क्रमः == CblasColMajor) {
		यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo = 0;
		यदि (Uplo == Cblas_अधस्तनम्) uplo = 1;

		यदि (TransA == CblasNoTrans)
			transa = 0;
		यदि (TransA == CblasTrans)
			transa = 1;
#यदि_अपरिभाषितम् COMPLEX
		यदि (TransA == CblasConjNoTrans)
			transa = 0;
		यदि (TransA == CblasConjTrans)
			transa = 1;
#अन्यथा
		यदि (TransA == CblasConjNoTrans)
			transa = 2;
		यदि (TransA == CblasConjTrans)
			transa = 3;
#यदि_अन्तः
		यदि (TransB == CblasNoTrans)
			transb = 0;
		यदि (TransB == CblasTrans)
			transb = 1;
#यदि_अपरिभाषितम् COMPLEX
		यदि (TransB == CblasConjNoTrans)
			transb = 0;
		यदि (TransB == CblasConjTrans)
			transb = 1;
#अन्यथा
		यदि (TransB == CblasConjNoTrans)
			transb = 2;
		यदि (TransB == CblasConjTrans)
			transb = 3;
#यदि_अन्तः

		a = (शून्य *)A;
		b = (शून्य *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		यदि (transa & 1) nrowa = k;
		nrowb = k;
		यदि (transb & 1) nrowb = m;

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

	यदि (क्रमः == CblasRowMajor) {

		a = (शून्य *)B;
		b = (शून्य *)A;

		lda = LDB;
		ldb = LDA;

		यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo = 1;
		यदि (Uplo == Cblas_अधस्तनम्) uplo = 0;

		यदि (TransB == CblasNoTrans)
			transa = 0;
		यदि (TransB == CblasTrans)
			transa = 1;
#यदि_अपरिभाषितम् COMPLEX
		यदि (TransB == CblasConjNoTrans)
			transa = 0;
		यदि (TransB == CblasConjTrans)
			transa = 1;
#अन्यथा
		यदि (TransB == CblasConjNoTrans)
			transa = 2;
		यदि (TransB == CblasConjTrans)
			transa = 3;
#यदि_अन्तः
		यदि (TransA == CblasNoTrans)
			transb = 0;
		यदि (TransA == CblasTrans)
			transb = 1;
#यदि_अपरिभाषितम् COMPLEX
		यदि (TransA == CblasConjNoTrans)
			transb = 0;
		यदि (TransA == CblasConjTrans)
			transb = 1;
#अन्यथा
		यदि (TransA == CblasConjNoTrans)
			transb = 2;
		यदि (TransA == CblasConjTrans)
			transb = 3;
#यदि_अन्तः

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		यदि (transa & 1) ncola = k;
		ncolb = k;
		यदि (transb & 1) ncolb = m;

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
#यदि_अन्तः
#यदि_अपरिभाषितम् COMPLEX
#परिभाषय GEMMT_सारणी_सूचकांक ((uplo << 2) | (transb << 1) | transa)
#अन्यथा
#परिभाषय GEMMT_सारणी_सूचकांक ((uplo << 4) | (transb << 2) | transa)
#यदि_अन्तः

	blas_arg_t तर्काः;
	FLOAT *बफर, *sa, *sb;

#यदि_परिभाषितम् SMP
	द्विगुण MMK;
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
	पूर्णाङ्क बहुलकम् = BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
	पूर्णाङ्क बहुलकम् = BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
	पूर्णाङ्क बहुलकम् = BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
	पूर्णाङ्क बहुलकम् = BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
	पूर्णाङ्क बहुलकम् = BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
	पूर्णाङ्क बहुलकम् = BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	स्थिर पूर्णाङ्क (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#यदि_अपरिभाषितम् COMPLEX
		GEMMT_UNN, GEMMT_UTN,
		GEMMT_UNT, GEMMT_UTT,
		GEMMT_LNN, GEMMT_LTN,
		GEMMT_LNT, GEMMT_LTT,
#अन्यथा
		GEMMT_UNN, GEMMT_UTN, GEMMT_URN, GEMMT_UCN,
		GEMMT_UNT, GEMMT_UTT, GEMMT_URT, GEMMT_UCT,
		GEMMT_UNR, GEMMT_UTR, GEMMT_URR, GEMMT_UCR,
		GEMMT_UNC, GEMMT_UTC, GEMMT_URC, GEMMT_UCC,
		GEMMT_LNN, GEMMT_LTN, GEMMT_LRN, GEMMT_LCN,
		GEMMT_LNT, GEMMT_LTT, GEMMT_LRT, GEMMT_LCT,
		GEMMT_LNR, GEMMT_LTR, GEMMT_LRR, GEMMT_LCR,
		GEMMT_LNC, GEMMT_LTC, GEMMT_LRC, GEMMT_LCC,
#यदि_अन्तः
	};

	यदि (m == 0)
		प्रत्यागम;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	तर्काः.m = m;
	तर्काः.n = m;
	तर्काः.k = k;

	तर्काः.a = (शून्य *)a;
	तर्काः.b = (शून्य *)b;
	तर्काः.c = (शून्य *)c;

	तर्काः.lda = lda;
	तर्काः.ldb = ldb;
	तर्काः.ldc = ldc;

#यदि defined(COMPLEX)
	तर्काः.alpha = (शून्य *)alpha;
	तर्काः.beta  = (शून्य *)beta;
#अन्यथा
	तर्काः.alpha = (शून्य *)&alpha;
	तर्काः.beta  = (शून्य *)&beta;
#यदि_अन्तः

#यदि_परिभाषितम् SMP
	तर्काः.common   = रिक्त;
	तर्काः.nthreads = 1;
#यदि_अन्तः

	बफर = (FLOAT *)blas_स्मृतिः_alloc(0);
	यदि (!बफर) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, परिमाण(ERROR_NAME));
		प्रत्यागम;
	}

	sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);

#यदि_परिभाषितम् SMP
	MMK = (द्विगुण)(m + 1) * (द्विगुण)m * (द्विगुण)k;
	यदि (MMK <= (SMP_देहली_MIN * GEMM_MULTITHREAD_देहली)) {
		तर्काः.nthreads = 1;
	} अन्यथा {
		तर्काः.nthreads = अङ्कः_cpu_avail(3);
	}

	यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

		(gemmt[GEMMT_सारणी_सूचकांक])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP
	} अन्यथा {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		बहुलकम् |= (uplo << BLAS_UPLO_सरणम्);

		syrk_thread(बहुलकम्, &तर्काः, रिक्त, रिक्त, gemmt[GEMMT_सारणी_सूचकांक], sa, sb, तर्काः.nthreads);
	}
#यदि_अन्तः

	blas_स्मृतिः_free(बफर);

	FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, तर्काः.m * तर्काः.k + तर्काः.k * तर्काः.n + तर्काः.m * तर्काः.n / 2,
			     तर्काः.m * तर्काः.n * तर्काः.k);

	IDEBUG_अन्तः;

	प्रत्यागम;
}
