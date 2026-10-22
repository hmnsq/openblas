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
#define ERROR_NAME "QGBMV "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#وگرنه
#define ERROR_NAME "SGBMV "
#endif

ایستا تهی (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, تهی *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#وگرنهاگر defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#وگرنه
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
ایستا صحیح (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#وگرنهاگر defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#وگرنه
      sgbmv_thread_n, sgbmv_thread_t,
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
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order,
	   شمارش CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  صحیح trans;
#ifdef SMP
  صحیح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  اگر (order == CblasColMajor) {
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

  اگر (order == CblasRowMajor) {
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
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#endif

  اگر ((m==0) || (n==0)) بازگشت;

  lenx = n;
  leny = m;
  اگر (trans) lenx = m;
  اگر (trans) leny = n;

  اگر (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), پوچ, 0, پوچ, 0);

  اگر (alpha == ZERO) بازگشت;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  اگر (incx < 0) x -= (lenx-1)*incx;
  اگر (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  اگر (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#ifdef SMP
  اگر (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  وگرنه
    nthreads = رقم_cpu_avail(2);

  اگر (nthreads == 1) {
#endif

  (gbmv[(صحیح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } وگرنه {

    (gbmv_thread[(صحیح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  بازگشت;
}
