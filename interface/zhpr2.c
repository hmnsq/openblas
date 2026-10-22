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
#घोषणा ERROR_NAME "XHPR2 "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZHPR2 "
#अन्यथा
#घोषणा ERROR_NAME "CHPR2 "
#पूर्ण_अगर

अटल पूर्णांक (*hpr2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#अगर_घोषित XDOUBLE
  xhpr2_U, xhpr2_L, xhpr2_V, xhpr2_M,
#या_अगर defined(DOUBLE)
  zhpr2_U, zhpr2_L, zhpr2_V, zhpr2_M,
#अन्यथा
  chpr2_U, chpr2_L, chpr2_V, chpr2_M,
#पूर्ण_अगर
};

#अगर_घोषित SMP
अटल पूर्णांक (*hpr2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, पूर्णांक) = {
#अगर_घोषित XDOUBLE
  xhpr2_thread_U, xhpr2_thread_L, xhpr2_thread_V, xhpr2_thread_M,
#या_अगर defined(DOUBLE)
  zhpr2_thread_U, zhpr2_thread_L, zhpr2_thread_V, zhpr2_thread_M,
#अन्यथा
  chpr2_thread_U, chpr2_thread_L, chpr2_thread_V, chpr2_thread_M,
#पूर्ण_अगर
};
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  अक्षर uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  पूर्णांक uplo;
  FLOAT *बफर;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  अगर (uplo_arg  == 'U') uplo  = 0;
  अगर (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  अगर (incy == 0)          info =  7;
  अगर (incx == 0)          info =  5;
  अगर (n < 0)              info =  2;
  अगर (uplo  < 0)          info =  1;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   क्रमागत CBLAS_UPLO Uplo,
	   blasint n,
	   व्योम *VALPHA,
	   व्योम  *vx, blasint incx,
	   व्योम  *vy, blasint incy,
	   व्योम  *va) {

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;
  FLOAT* a = (FLOAT*) va;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT *बफर;
  पूर्णांक uplo;
  blasint info;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  अगर (order == CblasColMajor) {
    अगर (Uplo == CblasUpper)         uplo  = 0;
    अगर (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    अगर (incy == 0)          info =  7;
    अगर (incx == 0)          info =  5;
    अगर (n < 0)              info =  2;
    अगर (uplo  < 0)          info =  1;
  }

  अगर (order == CblasRowMajor) {
    अगर (Uplo == CblasUpper)         uplo  = 3;
    अगर (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    अगर (incx == 0)          info =  7;
    अगर (incy == 0)          info =  5;
    अगर (n < 0)              info =  2;
    अगर (uplo  < 0)          info =  1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर

  अगर (n == 0) वापस;

  अगर ((alpha_r == ZERO) && (alpha_i == ZERO)) वापस;


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0 ) x -= (n - 1) * incx * 2;
  अगर (incy < 0 ) y -= (n - 1) * incy * 2;

  बफर = (FLOAT *)blas_memory_alloc(1);
  अगर (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अगर_घोषित SMP
  अगर (n < 100)
	nthreads = 1;
  अन्यथा
    nthreads = अंक_cpu_avail(2);

  अगर (nthreads == 1) {
#पूर्ण_अगर

    (hpr2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, बफर);

#अगर_घोषित SMP
  } अन्यथा {

    (hpr2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, बफर, nthreads);

  }
#पूर्ण_अगर

  blas_memory_free(बफर);

  FUNCTION_PROFILE_END(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  वापस;
}
