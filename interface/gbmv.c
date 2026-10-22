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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QGBMV "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DGBMV "
#وإلا
#تعریف_کړه ERROR_NAME "SGBMV "
#د_که_پای

ساكن فراغ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#که_تعریف_شوی XDOUBLE
      qgbmv_n, qgbmv_t,
#وإلاإذا defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#وإلا
      sgbmv_n, sgbmv_t,
#د_که_پای
};

#که_تعریف_شوی SMP
ساكن صحيح (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#که_تعریف_شوی XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#وإلاإذا defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#وإلا
      sgbmv_thread_n, sgbmv_thread_t,
#د_که_پای
};
#د_که_پای

OPENBLAS_EXPORT
#که_نه_تعریف_شوی CBLAS

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
  FLOAT *بفر;
#که_تعریف_شوی SMP
  صحيح nthreads;
#د_که_پای

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  إذا (trans == 'N') i = 0;
  إذا (trans == 'T') i = 1;
  إذا (trans == 'R') i = 0;
  إذا (trans == 'C') i = 1;

  إذا (incy == 0)	 info = 13;
  إذا (incx == 0)	 info = 10;
  إذا (lda < kl + ku + 1) info = 8;
  إذا (kl < 0)		 info = 5;
  إذا (ku < 0)		 info = 4;
  إذا (n < 0)		 info = 3;
  إذا (m < 0)		 info = 2;
  إذا (i < 0)		 info = 1;

  trans = i;

  إذا (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب,
	   تعداد CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *بفر;
  blasint lenx, leny, info, t;
  صحيح trans;
#که_تعریف_شوی SMP
  صحيح nthreads;
#د_که_پای

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (TransA == CblasNoTrans)     trans = 0;
    إذا (TransA == CblasTrans)       trans = 1;
    إذا (TransA == CblasConjNoTrans) trans = 0;
    إذا (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    إذا (incy == 0)	 info = 13;
    إذا (incx == 0)	 info = 10;
    إذا (lda < kl + ku + 1) info = 8;
    إذا (kl < 0)		 info = 5;
    إذا (ku < 0)		 info = 4;
    إذا (n < 0)		 info = 3;
    إذا (m < 0)		 info = 2;
    إذا (trans < 0)	 info = 1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    إذا (TransA == CblasNoTrans)     trans = 1;
    إذا (TransA == CblasTrans)       trans = 0;
    إذا (TransA == CblasConjNoTrans) trans = 1;
    إذا (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    إذا (incy == 0)	 info = 13;
    إذا (incx == 0)	 info = 10;
    إذا (lda < kl + ku + 1) info = 8;
    إذا (kl < 0)		 info = 5;
    إذا (ku < 0)		 info = 4;
    إذا (n < 0)		 info = 3;
    إذا (m < 0)		 info = 2;
    إذا (trans < 0)	 info = 1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#د_که_پای

  إذا ((m==0) || (n==0)) أرجع;

  lenx = n;
  leny = m;
  إذا (trans) lenx = m;
  إذا (trans) leny = n;

  إذا (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);

  إذا (alpha == ZERO) أرجع;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0) x -= (lenx-1)*incx;
  إذا (incy < 0) y -= (leny-1)*incy;

  بفر = (FLOAT *)blas_حافظه_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#که_تعریف_شوی SMP
  إذا (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  وإلا
    nthreads = عدد_cpu_avail(2);

  إذا (nthreads == 1) {
#د_که_پای

  (gbmv[(صحيح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, بفر);

#که_تعریف_شوی SMP
  } وإلا {

    (gbmv_thread[(صحيح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, بفر, nthreads);

  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(1, m * n / 2 + n, m * n);

  IDEBUG_پای_ټکی;

  أرجع;
}
