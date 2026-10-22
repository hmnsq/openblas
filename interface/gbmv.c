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
#define ERROR_NAME "QGBMV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#אחרת
#define ERROR_NAME "SGBMV "
#endif

ساكن فراغ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#وإلاإذا defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#אחרת
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
ساكن صحيح (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#وإلاإذا defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#אחרת
      sgbmv_thread_n, sgbmv_thread_t,
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

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  אם (trans == 'N') i = 0;
  אם (trans == 'T') i = 1;
  אם (trans == 'R') i = 0;
  אם (trans == 'C') i = 1;

  אם (incy == 0)	 info = 13;
  אם (incx == 0)	 info = 10;
  אם (lda < kl + ku + 1) info = 8;
  אם (kl < 0)		 info = 5;
  אם (ku < 0)		 info = 4;
  אם (n < 0)		 info = 3;
  אם (m < 0)		 info = 2;
  אם (i < 0)		 info = 1;

  trans = i;

  אם (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#אחרת

فراغ CNAME(تعداد CBLAS_ORDER order,
	   تعداد CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  صحيح trans;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  אם (order == CblasColMajor) {
    אם (TransA == CblasNoTrans)     trans = 0;
    אם (TransA == CblasTrans)       trans = 1;
    אם (TransA == CblasConjNoTrans) trans = 0;
    אם (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    אם (incy == 0)	 info = 13;
    אם (incx == 0)	 info = 10;
    אם (lda < kl + ku + 1) info = 8;
    אם (kl < 0)		 info = 5;
    אם (ku < 0)		 info = 4;
    אם (n < 0)		 info = 3;
    אם (m < 0)		 info = 2;
    אם (trans < 0)	 info = 1;
  }

  אם (order == CblasRowMajor) {
    אם (TransA == CblasNoTrans)     trans = 1;
    אם (TransA == CblasTrans)       trans = 0;
    אם (TransA == CblasConjNoTrans) trans = 1;
    אם (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    אם (incy == 0)	 info = 13;
    אם (incx == 0)	 info = 10;
    אם (lda < kl + ku + 1) info = 8;
    אם (kl < 0)		 info = 5;
    אם (ku < 0)		 info = 4;
    אם (n < 0)		 info = 3;
    אם (m < 0)		 info = 2;
    אם (trans < 0)	 info = 1;
  }

  אם (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#endif

  אם ((m==0) || (n==0)) החזר;

  lenx = n;
  leny = m;
  אם (trans) lenx = m;
  אם (trans) leny = n;

  אם (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);

  אם (alpha == ZERO) החזר;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  אם (incx < 0) x -= (lenx-1)*incx;
  אם (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  אם (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#ifdef SMP
  אם (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  אחרת
    nthreads = num_cpu_avail(2);

  אם (nthreads == 1) {
#endif

  (gbmv[(صحيح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } אחרת {

    (gbmv_thread[(صحيح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  החזר;
}
