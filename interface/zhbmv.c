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
#समावेश <ctype.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XHBMV "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZHBMV "
#अन्यथा
#घोषणा ERROR_NAME "CHBMV "
#पूर्ण_अगर

अटल  पूर्णांक (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, व्योम *) = {
#अगर_घोषित XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#या_अगर defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#अन्यथा
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#पूर्ण_अगर
};

#अगर_घोषित SMPBUG
अटल  पूर्णांक (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#अगर_घोषित XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#या_अगर defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#अन्यथा
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#पूर्ण_अगर
};
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  अक्षर uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  पूर्णांक uplo;
  FLOAT *बफर;
#अगर_घोषित SMPBUG
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  अगर (uplo_arg  == 'U') uplo  = 0;
  अगर (uplo_arg  == 'L') uplo  = 1;
  अगर (uplo_arg  == 'V') uplo  = 2;
  अगर (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  अगर (incy == 0)          info = 11;
  अगर (incx == 0)          info =  8;
  अगर (lda  < k + 1)       info =  6;
  अगर (k < 0)              info =  3;
  अगर (n < 0)              info =  2;
  अगर (uplo  < 0)          info =  1;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   क्रमागत CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   व्योम *VALPHA,
	   व्योम  *va, blasint lda,
	   व्योम  *vx, blasint incx,
	   व्योम *VBETA,
	   व्योम  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  FLOAT *बफर;
  पूर्णांक uplo;
  blasint info;
#अगर_घोषित SMPBUG
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  अगर (order == CblasColMajor) {
    अगर (Uplo == CblasUpper)         uplo  = 0;
    अगर (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    अगर (incy == 0)          info = 11;
    अगर (incx == 0)          info =  8;
    अगर (lda  < k + 1)       info =  6;
    अगर (k < 0)              info =  3;
    अगर (n < 0)              info =  2;
    अगर (uplo  < 0)          info =  1;
  }

  अगर (order == CblasRowMajor) {
    अगर (Uplo == CblasUpper)         uplo  = 3;
    अगर (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    अगर (incy == 0)          info = 11;
    अगर (incx == 0)          info =  8;
    अगर (lda  < k + 1)       info =  6;
    अगर (k < 0)              info =  3;
    अगर (n < 0)              info =  2;
    अगर (uplo  < 0)          info =  1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर

  अगर (n == 0) वापस;

  अगर ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), शून्य, 0, शून्य, 0);

  अगर ((alpha_r == ZERO) && (alpha_i == ZERO)) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  अगर (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  बफर = (FLOAT *)blas_memory_alloc(1);
  अगर (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अगर_घोषित SMPBUG
  nthreads = अंक_cpu_avail(2);

  अगर (nthreads == 1) {
#पूर्ण_अगर

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, बफर);

#अगर_घोषित SMPBUG
  } अन्यथा {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, बफर, nthreads);

  }
#पूर्ण_अगर

  blas_memory_free(बफर);

  FUNCTION_PROFILE_END(4, n * k / 2 + n,  n * k);

  IDEBUG_END;

  वापस;
}
