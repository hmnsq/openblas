/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
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
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGBMV "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGBMV "
#अन्यथा
#घोषणा ERROR_NAME "SGBMV "
#पूर्ण_अगर

अटल व्योम (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, व्योम *) = {
#अगर_घोषित XDOUBLE
      qgbmv_n, qgbmv_t,
#या_अगर defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#अन्यथा
      sgbmv_n, sgbmv_t,
#पूर्ण_अगर
};

#अगर_घोषित SMP
अटल पूर्णांक (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#अगर_घोषित XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#या_अगर defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#अन्यथा
      sgbmv_thread_n, sgbmv_thread_t,
#पूर्ण_अगर
};
#पूर्ण_अगर

OPENBLAS_EXPORT
#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  अक्षर trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *बफर;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  अगर (trans == 'N') i = 0;
  अगर (trans == 'T') i = 1;
  अगर (trans == 'R') i = 0;
  अगर (trans == 'C') i = 1;

  अगर (incy == 0)	 info = 13;
  अगर (incx == 0)	 info = 10;
  अगर (lda < kl + ku + 1) info = 8;
  अगर (kl < 0)		 info = 5;
  अगर (ku < 0)		 info = 4;
  अगर (n < 0)		 info = 3;
  अगर (m < 0)		 info = 2;
  अगर (i < 0)		 info = 1;

  trans = i;

  अगर (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   क्रमागत CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *बफर;
  blasint lenx, leny, info, t;
  पूर्णांक trans;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  अगर (order == CblasColMajor) {
    अगर (TransA == CblasNoTrans)     trans = 0;
    अगर (TransA == CblasTrans)       trans = 1;
    अगर (TransA == CblasConjNoTrans) trans = 0;
    अगर (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    अगर (incy == 0)	 info = 13;
    अगर (incx == 0)	 info = 10;
    अगर (lda < kl + ku + 1) info = 8;
    अगर (kl < 0)		 info = 5;
    अगर (ku < 0)		 info = 4;
    अगर (n < 0)		 info = 3;
    अगर (m < 0)		 info = 2;
    अगर (trans < 0)	 info = 1;
  }

  अगर (order == CblasRowMajor) {
    अगर (TransA == CblasNoTrans)     trans = 1;
    अगर (TransA == CblasTrans)       trans = 0;
    अगर (TransA == CblasConjNoTrans) trans = 1;
    अगर (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    अगर (incy == 0)	 info = 13;
    अगर (incx == 0)	 info = 10;
    अगर (lda < kl + ku + 1) info = 8;
    अगर (kl < 0)		 info = 5;
    अगर (ku < 0)		 info = 4;
    अगर (n < 0)		 info = 3;
    अगर (m < 0)		 info = 2;
    अगर (trans < 0)	 info = 1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर

  अगर ((m==0) || (n==0)) वापस;

  lenx = n;
  leny = m;
  अगर (trans) lenx = m;
  अगर (trans) leny = n;

  अगर (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), शून्य, 0, शून्य, 0);

  अगर (alpha == ZERO) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (lenx-1)*incx;
  अगर (incy < 0) y -= (leny-1)*incy;

  बफर = (FLOAT *)blas_memory_alloc(1);
  अगर (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अगर_घोषित SMP
  अगर (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  अन्यथा
    nthreads = अंक_cpu_avail(2);

  अगर (nthreads == 1) {
#पूर्ण_अगर

  (gbmv[(पूर्णांक)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, बफर);

#अगर_घोषित SMP
  } अन्यथा {

    (gbmv_thread[(पूर्णांक)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, बफर, nthreads);

  }
#पूर्ण_अगर

  blas_memory_free(बफर);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  वापस;
}
