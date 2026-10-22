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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGBMV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZGBMV "
#অন্যথা
#define ERROR_NAME "CGBMV "
#endif

স্থির শূন্য (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
      xgbmv_n, xgbmv_t, xgbmv_r, xgbmv_c,
      xgbmv_o, xgbmv_u, xgbmv_s, xgbmv_d,
#নোচেদ্যদি defined(DOUBLE)
      zgbmv_n, zgbmv_t, zgbmv_r, zgbmv_c,
      zgbmv_o, zgbmv_u, zgbmv_s, zgbmv_d,
#অন্যথা
      cgbmv_n, cgbmv_t, cgbmv_r, cgbmv_c,
      cgbmv_o, cgbmv_u, cgbmv_s, cgbmv_d,
#endif
};

#ifdef SMP
স্থির পূর্ণাঙ্ক (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT *,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
      xgbmv_thread_n, xgbmv_thread_t, xgbmv_thread_r, xgbmv_thread_c,
      xgbmv_thread_o, xgbmv_thread_u, xgbmv_thread_s, xgbmv_thread_d,
#নোচেদ্যদি defined(DOUBLE)
      zgbmv_thread_n, zgbmv_thread_t, zgbmv_thread_r, zgbmv_thread_c,
      zgbmv_thread_o, zgbmv_thread_u, zgbmv_thread_s, zgbmv_thread_d,
#অন্যথা
      cgbmv_thread_n, cgbmv_thread_t, cgbmv_thread_r, cgbmv_thread_c,
      cgbmv_thread_o, cgbmv_thread_u, cgbmv_thread_s, cgbmv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(ৱর্ণ *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  ৱর্ণ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
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

  যদি (trans == 'N')  i = 0;
  যদি (trans == 'T')  i = 1;
  যদি (trans == 'R')  i = 2;
  যদি (trans == 'C')  i = 3;
  যদি (trans == 'O')  i = 4;
  যদি (trans == 'U')  i = 5;
  যদি (trans == 'S')  i = 6;
  যদি (trans == 'D')  i = 7;

  যদি (incy == 0)	 info = 13;
  যদি (incx == 0)	 info = 10;
  যদি (lda < kl + ku + 1) info = 8;
  যদি (kl < 0)		 info = 5;
  যদি (ku < 0)		 info = 4;
  যদি (n < 0)		 info = 3;
  যদি (m < 0)		 info = 2;
  যদি (i < 0)		 info = 1;

  trans = i;

  যদি (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order,
	   গণনা CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   শূন্য *VALPHA,
	   শূন্য  *va, blasint lda,
	   শূন্য  *vx, blasint incx,
	   শূন্য *VBETA,
	   শূন্য  *vy, blasint incy){

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
  পূর্ণাঙ্ক trans;
  blasint info, t;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (TransA == CblasNoTrans)     trans = 0;
    যদি (TransA == CblasTrans)       trans = 1;
    যদি (TransA == CblasConjNoTrans) trans = 2;
    যদি (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    যদি (incy == 0)	 info = 13;
    যদি (incx == 0)	 info = 10;
    যদি (lda < kl + ku + 1) info = 8;
    যদি (kl < 0)		 info = 5;
    যদি (ku < 0)		 info = 4;
    যদি (n < 0)		 info = 3;
    যদি (m < 0)		 info = 2;
    যদি (trans < 0)	 info = 1;
  }

  যদি (order == CblasRowMajor) {
    যদি (TransA == CblasNoTrans)     trans = 1;
    যদি (TransA == CblasTrans)       trans = 0;
    যদি (TransA == CblasConjNoTrans) trans = 3;
    যদি (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    যদি (incy == 0)	 info = 13;
    যদি (incx == 0)	 info = 10;
    যদি (lda < kl + ku + 1) info = 8;
    যদি (kl < 0)		 info = 5;
    যদি (ku < 0)		 info = 4;
    যদি (n < 0)		 info = 3;
    যদি (m < 0)		 info = 2;
    যদি (trans < 0)	 info = 1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }

#endif

  যদি ((m==0) || (n==0)) ঘূৰাই;

  lenx = n;
  leny = m;
  যদি (trans & 1) lenx = m;
  যদি (trans & 1) leny = n;

  যদি (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), রিক্ত, 0, রিক্ত, 0);

  যদি (alpha_r == ZERO && alpha_i == ZERO) ঘূৰাই;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (lenx - 1) * incx * 2;
  যদি (incy < 0) y -= (leny - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }

#ifdef SMP
  যদি (m * n  < 125000 || ku + kl < 15) 
    nthreads = 1;
  অন্যথা
    nthreads = অঙ্ক_cpu_avail(2);

  যদি (nthreads == 1) {
#endif

  (gbmv[(পূর্ণাঙ্ক)trans])(m, n, kl, ku, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP

  } অন্যথা {

    (gbmv_thread[(পূর্ণাঙ্ক)trans])(m, n, kl, ku, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, m * n / 2 + n, m * n);

  IDEBUG_END;

  ঘূৰাই;
}
