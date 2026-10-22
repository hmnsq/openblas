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

#ضمن <stdio.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGBMV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZGBMV "
#ܐܚܪܢܐ
#define ERROR_NAME "CGBMV "
#endif

ساكن فراغ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
      xgbmv_n, xgbmv_t, xgbmv_r, xgbmv_c,
      xgbmv_o, xgbmv_u, xgbmv_s, xgbmv_d,
#وإلاإذا defined(DOUBLE)
      zgbmv_n, zgbmv_t, zgbmv_r, zgbmv_c,
      zgbmv_o, zgbmv_u, zgbmv_s, zgbmv_d,
#ܐܚܪܢܐ
      cgbmv_n, cgbmv_t, cgbmv_r, cgbmv_c,
      cgbmv_o, cgbmv_u, cgbmv_s, cgbmv_d,
#endif
};

#ifdef SMP
ساكن صحيح (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT *,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#ifdef XDOUBLE
      xgbmv_thread_n, xgbmv_thread_t, xgbmv_thread_r, xgbmv_thread_c,
      xgbmv_thread_o, xgbmv_thread_u, xgbmv_thread_s, xgbmv_thread_d,
#وإلاإذا defined(DOUBLE)
      zgbmv_thread_n, zgbmv_thread_t, zgbmv_thread_r, zgbmv_thread_c,
      zgbmv_thread_o, zgbmv_thread_u, zgbmv_thread_s, zgbmv_thread_d,
#ܐܚܪܢܐ
      cgbmv_thread_n, cgbmv_thread_t, cgbmv_thread_r, cgbmv_thread_c,
      cgbmv_thread_o, cgbmv_thread_u, cgbmv_thread_s, cgbmv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  حرف trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  صحيح nthreads;
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

  ܐܢ (trans == 'N')  i = 0;
  ܐܢ (trans == 'T')  i = 1;
  ܐܢ (trans == 'R')  i = 2;
  ܐܢ (trans == 'C')  i = 3;
  ܐܢ (trans == 'O')  i = 4;
  ܐܢ (trans == 'U')  i = 5;
  ܐܢ (trans == 'S')  i = 6;
  ܐܢ (trans == 'D')  i = 7;

  ܐܢ (incy == 0)	 info = 13;
  ܐܢ (incx == 0)	 info = 10;
  ܐܢ (lda < kl + ku + 1) info = 8;
  ܐܢ (kl < 0)		 info = 5;
  ܐܢ (ku < 0)		 info = 4;
  ܐܢ (n < 0)		 info = 3;
  ܐܢ (m < 0)		 info = 2;
  ܐܢ (i < 0)		 info = 1;

  trans = i;

  ܐܢ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ܐܚܪܢܐ

فراغ CNAME(تعداد CBLAS_ORDER order,
	   تعداد CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   فراغ *VALPHA,
	   فراغ  *va, blasint lda,
	   فراغ  *vx, blasint incx,
	   فراغ *VBETA,
	   فراغ  *vy, blasint incy){

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
  صحيح trans;
  blasint info, t;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ܐܢ (order == CblasColMajor) {
    ܐܢ (TransA == CblasNoTrans)     trans = 0;
    ܐܢ (TransA == CblasTrans)       trans = 1;
    ܐܢ (TransA == CblasConjNoTrans) trans = 2;
    ܐܢ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ܐܢ (incy == 0)	 info = 13;
    ܐܢ (incx == 0)	 info = 10;
    ܐܢ (lda < kl + ku + 1) info = 8;
    ܐܢ (kl < 0)		 info = 5;
    ܐܢ (ku < 0)		 info = 4;
    ܐܢ (n < 0)		 info = 3;
    ܐܢ (m < 0)		 info = 2;
    ܐܢ (trans < 0)	 info = 1;
  }

  ܐܢ (order == CblasRowMajor) {
    ܐܢ (TransA == CblasNoTrans)     trans = 1;
    ܐܢ (TransA == CblasTrans)       trans = 0;
    ܐܢ (TransA == CblasConjNoTrans) trans = 3;
    ܐܢ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ܐܢ (incy == 0)	 info = 13;
    ܐܢ (incx == 0)	 info = 10;
    ܐܢ (lda < kl + ku + 1) info = 8;
    ܐܢ (kl < 0)		 info = 5;
    ܐܢ (ku < 0)		 info = 4;
    ܐܢ (n < 0)		 info = 3;
    ܐܢ (m < 0)		 info = 2;
    ܐܢ (trans < 0)	 info = 1;
  }

  ܐܢ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#endif

  ܐܢ ((m==0) || (n==0)) ܦܢܐ;

  lenx = n;
  leny = m;
  ܐܢ (trans & 1) lenx = m;
  ܐܢ (trans & 1) leny = n;

  ܐܢ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), عدم, 0, عدم, 0);

  ܐܢ (alpha_r == ZERO && alpha_i == ZERO) ܦܢܐ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ܐܢ (incx < 0) x -= (lenx - 1) * incx * 2;
  ܐܢ (incy < 0) y -= (leny - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ܐܢ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ifdef SMP
  ܐܢ (m * n  < 125000 || ku + kl < 15) 
    nthreads = 1;
  ܐܚܪܢܐ
    nthreads = num_cpu_avail(2);

  ܐܢ (nthreads == 1) {
#endif

  (gbmv[(صحيح)trans])(m, n, kl, ku, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP

  } ܐܚܪܢܐ {

    (gbmv_thread[(صحيح)trans])(m, n, kl, ku, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, m * n / 2 + n, m * n);

  IDEBUG_END;

  ܦܢܐ;
}
