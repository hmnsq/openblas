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

#شامل <stdio.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGBMV "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZGBMV "
#وگرنه
#define ERROR_NAME "CGBMV "
#endif

ایستا تهی (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, تهی *) = {
#ifdef XDOUBLE
      xgbmv_n, xgbmv_t, xgbmv_r, xgbmv_c,
      xgbmv_o, xgbmv_u, xgbmv_s, xgbmv_d,
#وگرنهاگر defined(DOUBLE)
      zgbmv_n, zgbmv_t, zgbmv_r, zgbmv_c,
      zgbmv_o, zgbmv_u, zgbmv_s, zgbmv_d,
#وگرنه
      cgbmv_n, cgbmv_t, cgbmv_r, cgbmv_c,
      cgbmv_o, cgbmv_u, cgbmv_s, cgbmv_d,
#endif
};

#ifdef SMP
ایستا صحیح (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT *,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#ifdef XDOUBLE
      xgbmv_thread_n, xgbmv_thread_t, xgbmv_thread_r, xgbmv_thread_c,
      xgbmv_thread_o, xgbmv_thread_u, xgbmv_thread_s, xgbmv_thread_d,
#وگرنهاگر defined(DOUBLE)
      zgbmv_thread_n, zgbmv_thread_t, zgbmv_thread_r, zgbmv_thread_c,
      zgbmv_thread_o, zgbmv_thread_u, zgbmv_thread_s, zgbmv_thread_d,
#وگرنه
      cgbmv_thread_n, cgbmv_thread_t, cgbmv_thread_r, cgbmv_thread_c,
      cgbmv_thread_o, cgbmv_thread_u, cgbmv_thread_s, cgbmv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  نویسه trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  صحیح nthreads;
#endif

  FLOAT alpha_r = ALPHA[0];
  FLOAT alpha_i = ALPHA[1];
  FLOAT beta_r  = BETA[0];
  FLOAT beta_i  = BETA[1];

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  اگر (trans == 'N')  i = 0;
  اگر (trans == 'T')  i = 1;
  اگر (trans == 'R')  i = 2;
  اگر (trans == 'C')  i = 3;
  اگر (trans == 'O')  i = 4;
  اگر (trans == 'U')  i = 5;
  اگر (trans == 'S')  i = 6;
  اگر (trans == 'D')  i = 7;

  اگر (incy == 0)	 info = 13;
  اگر (incx == 0)	 info = 10;
  اگر (lda < kl + ku + 1) info = 8;
  اگر (kl < 0)		 info = 5;
  اگر (ku < 0)		 info = 4;
  اگر (n < 0)		 info = 3;
  اگر (m < 0)		 info = 2;
  اگر (i < 0)		 info = 1;

  trans = i;

  اگر (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order,
	   شمارش CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   تهی *VALPHA,
	   تهی  *va, blasint lda,
	   تهی  *vx, blasint incx,
	   تهی *VBETA,
	   تهی  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r = ALPHA[0];
  FLOAT alpha_i = ALPHA[1];
  FLOAT beta_r  = BETA[0];
  FLOAT beta_i  = BETA[1];

  FLOAT *buffer;
  blasint lenx, leny;
  صحیح trans;
  blasint info, t;
#ifdef SMP
  صحیح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  اگر (order == CblasColMajor) {
    اگر (TransA == CblasNoTrans)     trans = 0;
    اگر (TransA == CblasTrans)       trans = 1;
    اگر (TransA == CblasConjNoTrans) trans = 2;
    اگر (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    اگر (incy == 0)	 info = 13;
    اگر (incx == 0)	 info = 10;
    اگر (lda < kl + ku + 1) info = 8;
    اگر (kl < 0)		 info = 5;
    اگر (ku < 0)		 info = 4;
    اگر (n < 0)		 info = 3;
    اگر (m < 0)		 info = 2;
    اگر (trans < 0)	 info = 1;
  }

  اگر (order == CblasRowMajor) {
    اگر (TransA == CblasNoTrans)     trans = 1;
    اگر (TransA == CblasTrans)       trans = 0;
    اگر (TransA == CblasConjNoTrans) trans = 3;
    اگر (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    اگر (incy == 0)	 info = 13;
    اگر (incx == 0)	 info = 10;
    اگر (lda < kl + ku + 1) info = 8;
    اگر (kl < 0)		 info = 5;
    اگر (ku < 0)		 info = 4;
    اگر (n < 0)		 info = 3;
    اگر (m < 0)		 info = 2;
    اگر (trans < 0)	 info = 1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#endif

  اگر ((m==0) || (n==0)) بازگشت;

  lenx = n;
  leny = m;
  اگر (trans & 1) lenx = m;
  اگر (trans & 1) leny = n;

  اگر (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), پوچ, 0, پوچ, 0);

  اگر (alpha_r == ZERO && alpha_i == ZERO) بازگشت;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  اگر (incx < 0) x -= (lenx - 1) * incx * 2;
  اگر (incy < 0) y -= (leny - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  اگر (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#ifdef SMP
  اگر (m * n  < 125000 || ku + kl < 15) 
    nthreads = 1;
  وگرنه
    nthreads = رقم_cpu_avail(2);

  اگر (nthreads == 1) {
#endif

  (gbmv[(صحیح)trans])(m, n, kl, ku, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP

  } وگرنه {

    (gbmv_thread[(صحیح)trans])(m, n, kl, ku, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, m * n / 2 + n, m * n);

  IDEBUG_END;

  بازگشت;
}
