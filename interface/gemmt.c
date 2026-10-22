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

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश "common.h"

#अगर_अघोषित COMPLEX
#घोषणा SMP_THRESHOLD_MIN 65536.0
#अगर_घोषित RNAME
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGEMMTR"
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGEMMTR"
#या_अगर defined(BFLOAT16)
#घोषणा ERROR_NAME "SBGEMMTR"
#अन्यथा
#घोषणा ERROR_NAME "SGEMMTR"
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGEMMT "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGEMMT "
#या_अगर defined(BFLOAT16)
#घोषणा ERROR_NAME "SBGEMMT "
#अन्यथा
#घोषणा ERROR_NAME "SGEMMT "
#पूर्ण_अगर
#पूर्ण_अगर
#अन्यथा
#घोषणा SMP_THRESHOLD_MIN 8192.0
#अगर_घोषित RNAME
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XGEMMTR"
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMMTR"
#अन्यथा
#घोषणा ERROR_NAME "CGEMMTR"
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XGEMMT "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMMT "
#अन्यथा
#घोषणा ERROR_NAME "CGEMMT "
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर_अघोषित GEMM_MULTITHREAD_THRESHOLD
#घोषणा GEMM_MULTITHREAD_THRESHOLD 4
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *UPLO, अक्षर *TRANSA, अक्षर *TRANSB,
	  blasint * M, blasint * K,
	  FLOAT * Alpha,
	  IFLOAT * a, blasint * ldA,
	  IFLOAT * b, blasint * ldB, FLOAT * Beta, FLOAT * c, blasint * ldC)
{

	blasint m, k;
	blasint lda, ldb, ldc;
	पूर्णांक transa, transb, uplo;
	blasint info;

	अक्षर transA, transB, Uplo;
	blasint nrowa, nrowb;
#अगर !defined(COMPLEX)
	FLOAT alpha, beta;
#पूर्ण_अगर

	PRINT_DEBUG_NAME;

	m = *M;
	k = *K;

#अगर defined(COMPLEX)
	FLOAT *alpha = Alpha;
	FLOAT *beta = Beta;
#अन्यथा
	alpha = *Alpha;
	beta = *Beta;
#पूर्ण_अगर

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

	अगर (transA == 'N')
		transa = 0;
	अगर (transA == 'T')
		transa = 1;
#अगर_अघोषित COMPLEX
	अगर (transA == 'R')
		transa = 0;
	अगर (transA == 'C')
		transa = 1;
#अन्यथा
	अगर (transA == 'R')
		transa = 2;
	अगर (transA == 'C')
		transa = 3;
#पूर्ण_अगर

	अगर (transB == 'N')
		transb = 0;
	अगर (transB == 'T')
		transb = 1;
#अगर_अघोषित COMPLEX
	अगर (transB == 'R')
		transb = 0;
	अगर (transB == 'C')
		transb = 1;
#अन्यथा
	अगर (transB == 'R')
		transb = 2;
	अगर (transB == 'C')
		transb = 3;
#पूर्ण_अगर

	अगर (Uplo == 'U')
		uplo = 0;
	अगर (Uplo == 'L')
		uplo = 1;
	
	nrowa = m;
	अगर (transa & 1) nrowa = k;
	nrowb = k;
	अगर (transb & 1) nrowb = m;

	info = 0;

	अगर (ldc < MAX(1, m))
		info = 13;
	अगर (ldb < MAX(1, nrowb))
		info = 10;
	अगर (lda < MAX(1, nrowa))
		info = 8;
	अगर (k < 0)
		info = 5;
	अगर (m < 0)
		info = 4;
	अगर (transb < 0)
		info = 3;
	अगर (transa < 0)
		info = 2;
	अगर (uplo < 0)
		info = 1;

	अगर (info != 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, माप(ERROR_NAME));
		वापस;
	}
#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_UPLO Uplo,
	   क्रमागत CBLAS_TRANSPOSE TransA, क्रमागत CBLAS_TRANSPOSE TransB, blasint m,
	   blasint k,
#अगर_अघोषित COMPLEX
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
#अन्यथा
	   व्योम *valpha,
	   व्योम *va, blasint LDA,
	   व्योम *vb, blasint LDB, व्योम *vbeta, व्योम *vc, blasint ldc)
{
	FLOAT *alpha = (FLOAT *) valpha;
	FLOAT *beta = (FLOAT *) vbeta;
	FLOAT *A = (FLOAT *) va;
	FLOAT *B = (FLOAT *) vb;
	FLOAT *c = (FLOAT *) vc;
#पूर्ण_अगर
	पूर्णांक transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	FLOAT *a, *b;

	PRINT_DEBUG_CNAME;

	uplo = -1;
	transa = -1;
	transb = -1;
	info = 0;

	अगर (order == CblasColMajor) {
		अगर (Uplo == CblasUpper) uplo = 0;
		अगर (Uplo == CblasLower) uplo = 1;

		अगर (TransA == CblasNoTrans)
			transa = 0;
		अगर (TransA == CblasTrans)
			transa = 1;
#अगर_अघोषित COMPLEX
		अगर (TransA == CblasConjNoTrans)
			transa = 0;
		अगर (TransA == CblasConjTrans)
			transa = 1;
#अन्यथा
		अगर (TransA == CblasConjNoTrans)
			transa = 2;
		अगर (TransA == CblasConjTrans)
			transa = 3;
#पूर्ण_अगर
		अगर (TransB == CblasNoTrans)
			transb = 0;
		अगर (TransB == CblasTrans)
			transb = 1;
#अगर_अघोषित COMPLEX
		अगर (TransB == CblasConjNoTrans)
			transb = 0;
		अगर (TransB == CblasConjTrans)
			transb = 1;
#अन्यथा
		अगर (TransB == CblasConjNoTrans)
			transb = 2;
		अगर (TransB == CblasConjTrans)
			transb = 3;
#पूर्ण_अगर

		a = (व्योम *)A;
		b = (व्योम *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa, nrowb;

		nrowa = m;
		अगर (transa & 1) nrowa = k;
		nrowb = k;
		अगर (transb & 1) nrowb = m;

		अगर (ldc < MAX(1, m))
			info = 13;
		अगर (ldb < MAX(1, nrowb))
			info = 10;
		अगर (lda < MAX(1, nrowa))
			info = 8;
		अगर (k < 0)
			info = 5;
		अगर (m < 0)
			info = 4;
		अगर (transb < 0)
			info = 3;
		अगर (transa < 0)
			info = 2;
		अगर (uplo < 0)
			info = 1;
	}

	अगर (order == CblasRowMajor) {

		a = (व्योम *)B;
		b = (व्योम *)A;

		lda = LDB;
		ldb = LDA;

		अगर (Uplo == CblasUpper) uplo = 1;
		अगर (Uplo == CblasLower) uplo = 0;

		अगर (TransB == CblasNoTrans)
			transa = 0;
		अगर (TransB == CblasTrans)
			transa = 1;
#अगर_अघोषित COMPLEX
		अगर (TransB == CblasConjNoTrans)
			transa = 0;
		अगर (TransB == CblasConjTrans)
			transa = 1;
#अन्यथा
		अगर (TransB == CblasConjNoTrans)
			transa = 2;
		अगर (TransB == CblasConjTrans)
			transa = 3;
#पूर्ण_अगर
		अगर (TransA == CblasNoTrans)
			transb = 0;
		अगर (TransA == CblasTrans)
			transb = 1;
#अगर_अघोषित COMPLEX
		अगर (TransA == CblasConjNoTrans)
			transb = 0;
		अगर (TransA == CblasConjTrans)
			transb = 1;
#अन्यथा
		अगर (TransA == CblasConjNoTrans)
			transb = 2;
		अगर (TransA == CblasConjTrans)
			transb = 3;
#पूर्ण_अगर

		info = -1;

		blasint ncola, ncolb;

		ncola = m;
		अगर (transa & 1) ncola = k;
		ncolb = k;
		अगर (transb & 1) ncolb = m;

		अगर (ldc < MAX(1,m))
			info = 13;
		अगर (ldb < MAX(1, ncolb))
			info = 8;
		अगर (lda < MAX(1, ncola))
			info = 10;
		अगर (k < 0)
			info = 5;
		अगर (m < 0)
			info = 4;
		अगर (transb < 0)
			info = 2;
		अगर (transa < 0)
			info = 3;
		अगर (uplo < 0)
			info = 1;
	}

	अगर (info >= 0) {
		BLASFUNC(xerbla) (ERROR_NAME, &info, माप(ERROR_NAME));
		वापस;
	}
#पूर्ण_अगर
#अगर_अघोषित COMPLEX
#घोषणा GEMMT_TABLE_सूचकांक ((uplo << 2) | (transb << 1) | transa)
#अन्यथा
#घोषणा GEMMT_TABLE_सूचकांक ((uplo << 4) | (transb << 2) | transa)
#पूर्ण_अगर

	blas_arg_t args;
	FLOAT *बफर, *sa, *sb;

#अगर_घोषित SMP
	द्विगुन MMK;
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
	पूर्णांक mode = BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
	पूर्णांक mode = BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
	पूर्णांक mode = BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
	पूर्णांक mode = BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
	पूर्णांक mode = BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
	पूर्णांक mode = BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

	/* [uplo][op(B)][op(A)], as in interface/gemm.c */
	अटल पूर्णांक (*gemmt[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#अगर_अघोषित COMPLEX
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
#पूर्ण_अगर
	};

	अगर (m == 0)
		वापस;

	IDEBUG_START;

	FUNCTION_PROFILE_START();

	args.m = m;
	args.n = m;
	args.k = k;

	args.a = (व्योम *)a;
	args.b = (व्योम *)b;
	args.c = (व्योम *)c;

	args.lda = lda;
	args.ldb = ldb;
	args.ldc = ldc;

#अगर defined(COMPLEX)
	args.alpha = (व्योम *)alpha;
	args.beta  = (व्योम *)beta;
#अन्यथा
	args.alpha = (व्योम *)&alpha;
	args.beta  = (व्योम *)&beta;
#पूर्ण_अगर

#अगर_घोषित SMP
	args.common   = शून्य;
	args.nthreads = 1;
#पूर्ण_अगर

	बफर = (FLOAT *)blas_memory_alloc(0);
	अगर (!बफर) {
		info = -999;
		BLASFUNC(xerbla) (ERROR_NAME, &info, माप(ERROR_NAME));
		वापस;
	}

	sa = (FLOAT *)((BLASLONG)बफर + GEMM_OFFSET_A);
	sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#अगर_घोषित SMP
	MMK = (द्विगुन)(m + 1) * (द्विगुन)m * (द्विगुन)k;
	अगर (MMK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
		args.nthreads = 1;
	} अन्यथा {
		args.nthreads = अंक_cpu_avail(3);
	}

	अगर (args.nthreads == 1) {
#पूर्ण_अगर

		(gemmt[GEMMT_TABLE_सूचकांक])(&args, शून्य, शून्य, sa, sb, 0);

#अगर_घोषित SMP
	} अन्यथा {

		/* Split the triangle into column ranges of equal area; every
		   thread then runs the single-threaded driver on its own range
		   and writes a disjoint set of columns of C. */
		mode |= (uplo << BLAS_UPLO_SHIFT);

		syrk_thread(mode, &args, शून्य, शून्य, gemmt[GEMMT_TABLE_सूचकांक], sa, sb, args.nthreads);
	}
#पूर्ण_अगर

	blas_memory_free(बफर);

	FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n / 2,
			     args.m * args.n * args.k);

	IDEBUG_END;

	वापस;
}
