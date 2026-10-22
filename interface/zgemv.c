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
#घोषणा ERROR_NAME "XGEMV "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMV "
#अन्यथा
#घोषणा ERROR_NAME "CGEMV "
#पूर्ण_अगर

#अगर_घोषित SMP
अटल पूर्णांक (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#अगर_घोषित XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#या_अगर defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#अन्यथा
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#पूर्ण_अगर
};
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  अक्षर trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *बफर;
  पूर्णांक बफर_आकार;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
		FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
		  GEMV_N, GEMV_T, GEMV_R, GEMV_C,
		  GEMV_O, GEMV_U, GEMV_S, GEMV_D,
		};

  blasint    info;
  blasint    lenx, leny;
  blasint    i;

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

  FLOAT beta_r  = *(BETA + 0);
  FLOAT beta_i  = *(BETA + 1);

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i    = -1;

  अगर (trans == 'N')  i = 0;
  अगर (trans == 'T')  i = 1;
  अगर (trans == 'R')  i = 2;
  अगर (trans == 'C')  i = 3;
  अगर (trans == 'O')  i = 4;
  अगर (trans == 'U')  i = 5;
  अगर (trans == 'S')  i = 6;
  अगर (trans == 'D')  i = 7;

  अगर (incy == 0)      info = 11;
  अगर (incx == 0)      info = 8;
  अगर (lda < MAX(1,m)) info = 6;
  अगर (n < 0) 	      info = 3;
  अगर (m < 0) 	      info = 2;
  अगर (i < 0)          info = 1;

  trans = i;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   क्रमागत CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   व्योम *VALPHA,
	   व्योम  *va, blasint lda,
	   व्योम  *vx, blasint incx,
	   व्योम *VBETA,
	   व्योम  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *बफर;
  blasint    lenx, leny;
  पूर्णांक trans, बफर_आकार;
  blasint info, t;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
	    FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
	      GEMV_N, GEMV_T, GEMV_R, GEMV_C,
	      GEMV_O, GEMV_U, GEMV_S, GEMV_D,
	    };

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

  FLOAT beta_r  = *(BETA + 0);
  FLOAT beta_i  = *(BETA + 1);

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  अगर (order == CblasColMajor) {
    अगर (TransA == CblasNoTrans)     trans = 0;
    अगर (TransA == CblasTrans)       trans = 1;
    अगर (TransA == CblasConjNoTrans) trans = 2;
    अगर (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    अगर (incy == 0)	  info = 11;
    अगर (incx == 0)	  info = 8;
    अगर (lda < MAX(1, m))  info = 6;
    अगर (n < 0)		  info = 3;
    अगर (m < 0)		  info = 2;
    अगर (trans < 0)        info = 1;

  }

  अगर (order == CblasRowMajor) {
    अगर (TransA == CblasNoTrans)     trans = 1;
    अगर (TransA == CblasTrans)       trans = 0;
    अगर (TransA == CblasConjNoTrans) trans = 3;
    अगर (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    अगर (incy == 0)	  info = 11;
    अगर (incx == 0)	  info = 8;
    अगर (lda < MAX(1, m))  info = 6;
    अगर (n < 0)		  info = 3;
    अगर (m < 0)		  info = 2;
    अगर (trans < 0)        info = 1;

  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर

  /*  Quick return if possible. */

  अगर (m == 0 || n == 0) वापस;

  lenx = n;
  leny = m;

  अगर (trans & 1) lenx = m;
  अगर (trans & 1) leny = n;

  अगर (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), शून्य, 0, शून्य, 0);

  अगर (alpha_r == ZERO && alpha_i == ZERO) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (lenx - 1) * incx * 2;
  अगर (incy < 0) y -= (leny - 1) * incy * 2;

  बफर_आकार = 2 * (m + n) + 128 / माप(FLOAT);
#अगर_घोषित WINDOWS_ABI
  बफर_आकार += 160 / माप(FLOAT) ;
#पूर्ण_अगर
  // for alignment
  बफर_आकार = (बफर_आकार + 3) & ~3;
  ढेर_ALLOC(बफर_आकार, FLOAT, बफर);

#अगर defined(ARCH_X86_64) && defined(MAX_ढेर_ALLOC) && MAX_ढेर_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  अगर(trans && ढेर_alloc_आकार)
    स_रखो(बफर, 0, MIN(बफर_आकार_2, माप(FLOAT) * बफर_आकार));
#पूर्ण_अगर

#अगर_घोषित SMP

#अगर defined(_WIN64) && defined(_M_ARM64)
  अगर (m*n > 25000000L) 
    nthreads = अंक_cpu_avail(4);
  अन्यथा
    nthreads = 1;
#अन्यथा
  अगर (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  अन्यथा
    nthreads = अंक_cpu_avail(2);
#पूर्ण_अगर  

  अगर (nthreads == 1) {
#पूर्ण_अगर  

    (gemv[(पूर्णांक)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, बफर);

#अगर_घोषित SMP
  } अन्यथा {
    (gemv_thread[(पूर्णांक)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, बफर, nthreads);
  }
#पूर्ण_अगर


  ढेर_FREE(बफर);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  वापस;
}
