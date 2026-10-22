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

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश "common.h"

#घोषणा SMP_THRESHOLD_MIN 65536.0
#घोषणा ERROR_NAME "SBGEMMT "

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
	IFLOAT *बफर;
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

	अगर (transA == 'N')
		transa = 0;
	अगर (transA == 'T')
		transa = 1;

	अगर (transA == 'R')
		transa = 0;
	अगर (transA == 'C')
		transa = 1;

	अगर (transB == 'N')
		transb = 0;
	अगर (transB == 'T')
		transb = 1;

	अगर (transB == 'R')
		transb = 0;
	अगर (transB == 'C')
		transb = 1;

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
	   FLOAT alpha,
	   IFLOAT * A, blasint LDA,
	   IFLOAT * B, blasint LDB, FLOAT beta, FLOAT * c, blasint ldc)
{
	IFLOAT *aa, *bb;
        FLOAT *cc;

	पूर्णांक transa, transb, uplo;
	blasint info;
	blasint lda, ldb;
	IFLOAT *a, *b;
	XFLOAT *बफर;

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

		अगर (TransA == CblasConjNoTrans)
			transa = 0;
		अगर (TransA == CblasConjTrans)
			transa = 1;

		अगर (TransB == CblasNoTrans)
			transb = 0;
		अगर (TransB == CblasTrans)
			transb = 1;

		अगर (TransB == CblasConjNoTrans)
			transb = 0;
		अगर (TransB == CblasConjTrans)
			transb = 1;

		a = (व्योम *)A;
		b = (व्योम *)B;
		lda = LDA;
		ldb = LDB;

		info = -1;

		blasint nrowa;
		blasint nrowb;
		nrowa = m;
		अगर (transa & 1) nrowa = k;
		nrowb = k;
		अगर (transb & 1)  nrowb = m;

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

		अगर (Uplo == CblasUpper) uplo = 0;
		अगर (Uplo == CblasLower) uplo = 1;

		अगर (TransB == CblasNoTrans)
			transa = 0;
		अगर (TransB == CblasTrans)
			transa = 1;

		अगर (TransB == CblasConjNoTrans)
			transa = 0;
		अगर (TransB == CblasConjTrans)
			transa = 1;

		अगर (TransA == CblasNoTrans)
			transb = 0;
		अगर (TransA == CblasTrans)
			transb = 1;

		अगर (TransA == CblasConjNoTrans)
			transb = 0;
		अगर (TransA == CblasConjTrans)
			transb = 1;

		info = -1;

		blasint ncola; 
		blasint ncolb;

		ncola = m;
		अगर (transa & 1) ncola = k;
		ncolb = k;

		अगर (transb & 1) {
			ncolb = m;
		}

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
	पूर्णांक बफर_आकार;
	blasint i, j;

#अगर_घोषित SMP
	पूर्णांक nthreads;
#पूर्ण_अगर


#अगर_घोषित SMP
	अटल पूर्णांक (*gemv_thread[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *,
				     BLASLONG, IFLOAT *, BLASLONG, FLOAT,
				     FLOAT *, BLASLONG, पूर्णांक) = {
		sbgemv_thread_n, sbgemv_thread_t,
	};
#पूर्ण_अगर
	पूर्णांक (*gemv[]) (BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,
		       IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
	GEMV_N, GEMV_T,};


	अगर (m == 0)
		वापस;

	IDEBUG_START;

	स्थिर blasint incb = ((transb & 1) == 0) ? 1 : ldb;

	अगर (uplo == 1) {
		क्रम (i = 0; i < m; i++) {
			j = m - i;

			aa = a + i;
			bb = b + i * ldb;
			अगर (transa & 1) {
				aa = a + lda * i;
			}
			अगर (transb & 1)
				bb = b + i;
			cc = c + i * ldc + i;

#अगर 0
			अगर (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, शून्य, 0, शून्य, 0);

			अगर (alpha == ZERO)
				जारी;
#पूर्ण_अगर

			IDEBUG_START;

			बफर_आकार = j + k + 128 / माप(FLOAT);
#अगर_घोषित WINDOWS_ABI
			बफर_आकार += 160 / माप(FLOAT);
#पूर्ण_अगर
			// for alignment
			बफर_आकार = (बफर_आकार + 3) & ~3;
			ढेर_ALLOC(बफर_आकार, IFLOAT, बफर);

#अगर_घोषित SMP

			अगर (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			अन्यथा
				nthreads = अंक_cpu_avail(2);

			अगर (nthreads == 1) {
#पूर्ण_अगर

				अगर (!(transa & 1))
				(gemv[(पूर्णांक)transa]) (j, k, alpha, aa, lda,
						     bb, incb, beta, cc, 1);
				अन्यथा
				(gemv[(पूर्णांक)transa]) (k, j, alpha, aa, lda,
						     bb, incb, beta, cc, 1);

#अगर_घोषित SMP
			} अन्यथा {
				अगर (!(transa & 1))
				(gemv_thread[(पूर्णांक)transa]) (j, k, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);
				अन्यथा
				(gemv_thread[(पूर्णांक)transa]) (k, j, alpha, aa,
							    lda, bb, incb, beta, cc,
							    1, nthreads);

			}
#पूर्ण_अगर

			ढेर_FREE(बफर);
		}
	} अन्यथा {

		क्रम (i = 0; i < m; i++) {
			j = i + 1;

			bb = b + i * ldb;
			अगर (transb & 1) {
				bb = b + i;
			}
			cc = c + i * ldc;

#अगर 0
			अगर (beta != ONE)
				SCAL_K(l, 0, 0, beta, cc, 1, शून्य, 0, शून्य, 0);

			अगर (alpha == ZERO)
				जारी;
#पूर्ण_अगर
			IDEBUG_START;

			बफर_आकार = j + k + 128 / माप(FLOAT);
#अगर_घोषित WINDOWS_ABI
			बफर_आकार += 160 / माप(FLOAT);
#पूर्ण_अगर
			// for alignment
			बफर_आकार = (बफर_आकार + 3) & ~3;
			ढेर_ALLOC(बफर_आकार, IFLOAT, बफर);

#अगर_घोषित SMP

			अगर (1L * j * k < 2304L * GEMM_MULTITHREAD_THRESHOLD)
				nthreads = 1;
			अन्यथा
				nthreads = अंक_cpu_avail(2);

			अगर (nthreads == 1) {
#पूर्ण_अगर

				अगर (!(transa & 1))
				(gemv[(पूर्णांक)transa]) (j, k, alpha, a, lda, bb,
						     incb, beta, cc, 1);
				अन्यथा
				(gemv[(पूर्णांक)transa]) (k, j, alpha, a, lda, bb,
						     incb, beta, cc, 1);

#अगर_घोषित SMP
			} अन्यथा {
				अगर (!(transa & 1))
				(gemv_thread[(पूर्णांक)transa]) (j, k, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
				अन्यथा
				(gemv_thread[(पूर्णांक)transa]) (k, j, alpha, a, lda,
							    bb, incb, beta, cc, 1,
							    nthreads);
			}
#पूर्ण_अगर

			ढेर_FREE(बफर);
		}
	}

	IDEBUG_END;

	वापस;
}
