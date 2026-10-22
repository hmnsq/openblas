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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGBMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGBMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGBMV "
#endif

ಸ್ಥಿರ ಶೂನ್ಯ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ಶೂನ್ಯ *) = {
#ifdef XDOUBLE
      xgbmv_n, xgbmv_t, xgbmv_r, xgbmv_c,
      xgbmv_o, xgbmv_u, xgbmv_s, xgbmv_d,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
      zgbmv_n, zgbmv_t, zgbmv_r, zgbmv_c,
      zgbmv_o, zgbmv_u, zgbmv_s, zgbmv_d,
#ಅನ್ಯಥಾ
      cgbmv_n, cgbmv_t, cgbmv_r, cgbmv_c,
      cgbmv_o, cgbmv_u, cgbmv_s, cgbmv_d,
#endif
};

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT *,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
      xgbmv_thread_n, xgbmv_thread_t, xgbmv_thread_r, xgbmv_thread_c,
      xgbmv_thread_o, xgbmv_thread_u, xgbmv_thread_s, xgbmv_thread_d,
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
      zgbmv_thread_n, zgbmv_thread_t, zgbmv_thread_r, zgbmv_thread_c,
      zgbmv_thread_o, zgbmv_thread_u, zgbmv_thread_s, zgbmv_thread_d,
#ಅನ್ಯಥಾ
      cgbmv_thread_n, cgbmv_thread_t, cgbmv_thread_r, cgbmv_thread_c,
      cgbmv_thread_o, cgbmv_thread_u, cgbmv_thread_s, cgbmv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  ವರ್ಣ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
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

  ಯದಿ (trans == 'N')  i = 0;
  ಯದಿ (trans == 'T')  i = 1;
  ಯದಿ (trans == 'R')  i = 2;
  ಯದಿ (trans == 'C')  i = 3;
  ಯದಿ (trans == 'O')  i = 4;
  ಯದಿ (trans == 'U')  i = 5;
  ಯದಿ (trans == 'S')  i = 6;
  ಯದಿ (trans == 'D')  i = 7;

  ಯದಿ (incy == 0)	 info = 13;
  ಯದಿ (incx == 0)	 info = 10;
  ಯದಿ (lda < kl + ku + 1) info = 8;
  ಯದಿ (kl < 0)		 info = 5;
  ಯದಿ (ku < 0)		 info = 4;
  ಯದಿ (n < 0)		 info = 3;
  ಯದಿ (m < 0)		 info = 2;
  ಯದಿ (i < 0)		 info = 1;

  trans = i;

  ಯದಿ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order,
	   ಗಣನಾ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   ಶೂನ್ಯ *VALPHA,
	   ಶೂನ್ಯ  *va, blasint lda,
	   ಶೂನ್ಯ  *vx, blasint incx,
	   ಶೂನ್ಯ *VBETA,
	   ಶೂನ್ಯ  *vy, blasint incy){

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
  ಪೂರ್ಣಾಙ್ಕ trans;
  blasint info, t;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (TransA == CblasNoTrans)     trans = 0;
    ಯದಿ (TransA == CblasTrans)       trans = 1;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 2;
    ಯದಿ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ಯದಿ (incy == 0)	 info = 13;
    ಯದಿ (incx == 0)	 info = 10;
    ಯದಿ (lda < kl + ku + 1) info = 8;
    ಯದಿ (kl < 0)		 info = 5;
    ಯದಿ (ku < 0)		 info = 4;
    ಯದಿ (n < 0)		 info = 3;
    ಯದಿ (m < 0)		 info = 2;
    ಯದಿ (trans < 0)	 info = 1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (TransA == CblasNoTrans)     trans = 1;
    ಯದಿ (TransA == CblasTrans)       trans = 0;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 3;
    ಯದಿ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ಯದಿ (incy == 0)	 info = 13;
    ಯದಿ (incx == 0)	 info = 10;
    ಯದಿ (lda < kl + ku + 1) info = 8;
    ಯದಿ (kl < 0)		 info = 5;
    ಯದಿ (ku < 0)		 info = 4;
    ಯದಿ (n < 0)		 info = 3;
    ಯದಿ (m < 0)		 info = 2;
    ಯದಿ (trans < 0)	 info = 1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  ಯದಿ ((m==0) || (n==0)) ಪ್ರತ್ಯಾಗಮ;

  lenx = n;
  leny = m;
  ಯದಿ (trans & 1) lenx = m;
  ಯದಿ (trans & 1) leny = n;

  ಯದಿ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ (alpha_r == ZERO && alpha_i == ZERO) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (lenx - 1) * incx * 2;
  ಯದಿ (incy < 0) y -= (leny - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ifdef SMP
  ಯದಿ (m * n  < 125000 || ku + kl < 15) 
    nthreads = 1;
  ಅನ್ಯಥಾ
    nthreads = ಅಂಕೆ_cpu_avail(2);

  ಯದಿ (nthreads == 1) {
#endif

  (gbmv[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, kl, ku, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP

  } ಅನ್ಯಥಾ {

    (gbmv_thread[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, kl, ku, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, m * n / 2 + n, m * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
