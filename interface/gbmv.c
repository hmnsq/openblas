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
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGBMV "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGBMV "
#ورنہ
#تعریف_کرو ERROR_NAME "SGBMV "
#اگر_ختم

ساکن خالی (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, خالی *) = {
#اگر_تعریف_شدہ XDOUBLE
      qgbmv_n, qgbmv_t,
#وگرنہ_اگر defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#ورنہ
      sgbmv_n, sgbmv_t,
#اگر_ختم
};

#اگر_تعریف_شدہ SMP
ساکن صحیح (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#اگر_تعریف_شدہ XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#وگرنہ_اگر defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#ورنہ
      sgbmv_thread_n, sgbmv_thread_t,
#اگر_ختم
};
#اگر_ختم

OPENBLAS_EXPORT
#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *TRANS, blasint *M, blasint *N,
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
  FLOAT *بفر;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  اگر (trans == 'N') i = 0;
  اگر (trans == 'T') i = 1;
  اگر (trans == 'R') i = 0;
  اگر (trans == 'C') i = 1;

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
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب,
	   شمار CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *بفر;
  blasint lenx, leny, info, t;
  صحیح trans;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    اگر (TransA == CblasNoTrans)     trans = 0;
    اگر (TransA == CblasTrans)       trans = 1;
    اگر (TransA == CblasConjNoTrans) trans = 0;
    اگر (TransA == CblasConjTrans)   trans = 1;

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

  اگر (ترتیب == Cblas_صف_Major) {
    اگر (TransA == CblasNoTrans)     trans = 1;
    اگر (TransA == CblasTrans)       trans = 0;
    اگر (TransA == CblasConjNoTrans) trans = 1;
    اگر (TransA == CblasConjTrans)   trans = 0;

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
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_ختم

  اگر ((m==0) || (n==0)) واپس;

  lenx = n;
  leny = m;
  اگر (trans) lenx = m;
  اگر (trans) leny = n;

  اگر (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), تہی, 0, تہی, 0);

  اگر (alpha == ZERO) واپس;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0) x -= (lenx-1)*incx;
  اگر (incy < 0) y -= (leny-1)*incy;

  بفر = (FLOAT *)blas_حافظہ_alloc(1);
  اگر (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_تعریف_شدہ SMP
  اگر (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  ورنہ
    nthreads = عدد_cpu_avail(2);

  اگر (nthreads == 1) {
#اگر_ختم

  (gbmv[(صحیح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, بفر);

#اگر_تعریف_شدہ SMP
  } ورنہ {

    (gbmv_thread[(صحیح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, بفر, nthreads);

  }
#اگر_ختم

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(1, m * n / 2 + n, m * n);

  IDEBUG_انتہا;

  واپس;
}
