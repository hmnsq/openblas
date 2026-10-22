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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#অন্যথা
#define ERROR_NAME "SGBMV "
#endif

স্থির শূন্য (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#নোচেদ্যদি defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#অন্যথা
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
স্থির পূর্ণাঙ্ক (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#নোচেদ্যদি defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#অন্যথা
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

শূন্য NAME(বর্ণ *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  বর্ণ trans = *TRANS;
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

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  যদি (trans == 'N') i = 0;
  যদি (trans == 'T') i = 1;
  যদি (trans == 'R') i = 0;
  যদি (trans == 'C') i = 1;

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
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order,
	   গণনা CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  পূর্ণাঙ্ক trans;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (TransA == CblasNoTrans)     trans = 0;
    যদি (TransA == CblasTrans)       trans = 1;
    যদি (TransA == CblasConjNoTrans) trans = 0;
    যদি (TransA == CblasConjTrans)   trans = 1;

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
    যদি (TransA == CblasConjNoTrans) trans = 1;
    যদি (TransA == CblasConjTrans)   trans = 0;

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
    প্রত্যাগম;
  }

#endif

  যদি ((m==0) || (n==0)) প্রত্যাগম;

  lenx = n;
  leny = m;
  যদি (trans) lenx = m;
  যদি (trans) leny = n;

  যদি (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), রিক্ত, 0, রিক্ত, 0);

  যদি (alpha == ZERO) প্রত্যাগম;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (lenx-1)*incx;
  যদি (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#ifdef SMP
  যদি (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  অন্যথা
    nthreads = অঙ্ক_cpu_avail(2);

  যদি (nthreads == 1) {
#endif

  (gbmv[(পূর্ণাঙ্ক)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } অন্যথা {

    (gbmv_thread[(পূর্ণাঙ্ক)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  প্রত্যাগম;
}
