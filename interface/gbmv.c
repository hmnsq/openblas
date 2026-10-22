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
#وإلا
#define ERROR_NAME "SGBMV "
#endif

ساكن فراغ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#وإلاإذا defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#وإلا
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
#وإلا
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

  𐤀𐤌 (trans == 'N') i = 0;
  𐤀𐤌 (trans == 'T') i = 1;
  𐤀𐤌 (trans == 'R') i = 0;
  𐤀𐤌 (trans == 'C') i = 1;

  𐤀𐤌 (incy == 0)	 info = 13;
  𐤀𐤌 (incx == 0)	 info = 10;
  𐤀𐤌 (lda < kl + ku + 1) info = 8;
  𐤀𐤌 (kl < 0)		 info = 5;
  𐤀𐤌 (ku < 0)		 info = 4;
  𐤀𐤌 (n < 0)		 info = 3;
  𐤀𐤌 (m < 0)		 info = 2;
  𐤀𐤌 (i < 0)		 info = 1;

  trans = i;

  𐤀𐤌 (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

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

  𐤀𐤌 (order == CblasColMajor) {
    𐤀𐤌 (TransA == CblasNoTrans)     trans = 0;
    𐤀𐤌 (TransA == CblasTrans)       trans = 1;
    𐤀𐤌 (TransA == CblasConjNoTrans) trans = 0;
    𐤀𐤌 (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    𐤀𐤌 (incy == 0)	 info = 13;
    𐤀𐤌 (incx == 0)	 info = 10;
    𐤀𐤌 (lda < kl + ku + 1) info = 8;
    𐤀𐤌 (kl < 0)		 info = 5;
    𐤀𐤌 (ku < 0)		 info = 4;
    𐤀𐤌 (n < 0)		 info = 3;
    𐤀𐤌 (m < 0)		 info = 2;
    𐤀𐤌 (trans < 0)	 info = 1;
  }

  𐤀𐤌 (order == CblasRowMajor) {
    𐤀𐤌 (TransA == CblasNoTrans)     trans = 1;
    𐤀𐤌 (TransA == CblasTrans)       trans = 0;
    𐤀𐤌 (TransA == CblasConjNoTrans) trans = 1;
    𐤀𐤌 (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    𐤀𐤌 (incy == 0)	 info = 13;
    𐤀𐤌 (incx == 0)	 info = 10;
    𐤀𐤌 (lda < kl + ku + 1) info = 8;
    𐤀𐤌 (kl < 0)		 info = 5;
    𐤀𐤌 (ku < 0)		 info = 4;
    𐤀𐤌 (n < 0)		 info = 3;
    𐤀𐤌 (m < 0)		 info = 2;
    𐤀𐤌 (trans < 0)	 info = 1;
  }

  𐤀𐤌 (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#endif

  𐤀𐤌 ((m==0) || (n==0)) أرجع;

  lenx = n;
  leny = m;
  𐤀𐤌 (trans) lenx = m;
  𐤀𐤌 (trans) leny = n;

  𐤀𐤌 (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);

  𐤀𐤌 (alpha == ZERO) أرجع;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  𐤀𐤌 (incx < 0) x -= (lenx-1)*incx;
  𐤀𐤌 (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  𐤀𐤌 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#ifdef SMP
  𐤀𐤌 (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  وإلا
    nthreads = num_cpu_avail(2);

  𐤀𐤌 (nthreads == 1) {
#endif

  (gbmv[(صحيح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } وإلا {

    (gbmv_thread[(صحيح)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  أرجع;
}
