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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGBMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGBMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGBMV "
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਖਾਲੀ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ਖਾਲੀ *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
      xgbmv_n, xgbmv_t, xgbmv_r, xgbmv_c,
      xgbmv_o, xgbmv_u, xgbmv_s, xgbmv_d,
#ਜਾਂਜੇ defined(DOUBLE)
      zgbmv_n, zgbmv_t, zgbmv_r, zgbmv_c,
      zgbmv_o, zgbmv_u, zgbmv_s, zgbmv_d,
#ਨਹੀਂਤਾਂ
      cgbmv_n, cgbmv_t, cgbmv_r, cgbmv_c,
      cgbmv_o, cgbmv_u, cgbmv_s, cgbmv_d,
#ਜੇ_ਅੰਤ
};

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT *,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
      xgbmv_thread_n, xgbmv_thread_t, xgbmv_thread_r, xgbmv_thread_c,
      xgbmv_thread_o, xgbmv_thread_u, xgbmv_thread_s, xgbmv_thread_d,
#ਜਾਂਜੇ defined(DOUBLE)
      zgbmv_thread_n, zgbmv_thread_t, zgbmv_thread_r, zgbmv_thread_c,
      zgbmv_thread_o, zgbmv_thread_u, zgbmv_thread_s, zgbmv_thread_d,
#ਨਹੀਂਤਾਂ
      cgbmv_thread_n, cgbmv_thread_t, cgbmv_thread_r, cgbmv_thread_c,
      cgbmv_thread_o, cgbmv_thread_u, cgbmv_thread_s, cgbmv_thread_d,
#ਜੇ_ਅੰਤ
};
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  ਅੱਖਰ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

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

  ਜੇ (trans == 'N')  i = 0;
  ਜੇ (trans == 'T')  i = 1;
  ਜੇ (trans == 'R')  i = 2;
  ਜੇ (trans == 'C')  i = 3;
  ਜੇ (trans == 'O')  i = 4;
  ਜੇ (trans == 'U')  i = 5;
  ਜੇ (trans == 'S')  i = 6;
  ਜੇ (trans == 'D')  i = 7;

  ਜੇ (incy == 0)	 info = 13;
  ਜੇ (incx == 0)	 info = 10;
  ਜੇ (lda < kl + ku + 1) info = 8;
  ਜੇ (kl < 0)		 info = 5;
  ਜੇ (ku < 0)		 info = 4;
  ਜੇ (n < 0)		 info = 3;
  ਜੇ (m < 0)		 info = 2;
  ਜੇ (i < 0)		 info = 1;

  trans = i;

  ਜੇ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ,
	   ਗਿਣਤੀ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   ਖਾਲੀ *VALPHA,
	   ਖਾਲੀ  *va, blasint lda,
	   ਖਾਲੀ  *vx, blasint incx,
	   ਖਾਲੀ *VBETA,
	   ਖਾਲੀ  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r = ALPHA[0];
  FLOAT alpha_i = ALPHA[1];
  FLOAT beta_r  = BETA[0];
  FLOAT beta_i  = BETA[1];

  FLOAT *ਬਫਰ;
  blasint lenx, leny;
  ਪੂਰਨ trans;
  blasint info, t;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਜੇ (TransA == CblasNoTrans)     trans = 0;
    ਜੇ (TransA == CblasTrans)       trans = 1;
    ਜੇ (TransA == CblasConjNoTrans) trans = 2;
    ਜੇ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ਜੇ (incy == 0)	 info = 13;
    ਜੇ (incx == 0)	 info = 10;
    ਜੇ (lda < kl + ku + 1) info = 8;
    ਜੇ (kl < 0)		 info = 5;
    ਜੇ (ku < 0)		 info = 4;
    ਜੇ (n < 0)		 info = 3;
    ਜੇ (m < 0)		 info = 2;
    ਜੇ (trans < 0)	 info = 1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਜੇ (TransA == CblasNoTrans)     trans = 1;
    ਜੇ (TransA == CblasTrans)       trans = 0;
    ਜੇ (TransA == CblasConjNoTrans) trans = 3;
    ਜੇ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ਜੇ (incy == 0)	 info = 13;
    ਜੇ (incx == 0)	 info = 10;
    ਜੇ (lda < kl + ku + 1) info = 8;
    ਜੇ (kl < 0)		 info = 5;
    ਜੇ (ku < 0)		 info = 4;
    ਜੇ (n < 0)		 info = 3;
    ਜੇ (m < 0)		 info = 2;
    ਜੇ (trans < 0)	 info = 1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  ਜੇ ((m==0) || (n==0)) ਵਾਪਸ;

  lenx = n;
  leny = m;
  ਜੇ (trans & 1) lenx = m;
  ਜੇ (trans & 1) leny = n;

  ਜੇ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ਨਿਰਾ, 0, ਨਿਰਾ, 0);

  ਜੇ (alpha_r == ZERO && alpha_i == ZERO) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0) x -= (lenx - 1) * incx * 2;
  ਜੇ (incy < 0) y -= (leny - 1) * incy * 2;

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਜੇ (m * n  < 125000 || ku + kl < 15) 
    nthreads = 1;
  ਨਹੀਂਤਾਂ
    nthreads = ਅੰਕ_cpu_avail(2);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  (gbmv[(ਪੂਰਨ)trans])(m, n, kl, ku, alpha_r, alpha_i, a, lda, x, incx, y, incy, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP

  } ਨਹੀਂਤਾਂ {

    (gbmv_thread[(ਪੂਰਨ)trans])(m, n, kl, ku, ALPHA, a, lda, x, incx, y, incy, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(4, m * n / 2 + n, m * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
