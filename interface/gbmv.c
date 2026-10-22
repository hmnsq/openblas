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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#ଅନ୍ୟଥା
#define ERROR_NAME "SGBMV "
#endif

ସ୍ଥିର ଶୂନ୍ୟ (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ଶୂନ୍ୟ *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#ନଚେତଯଦି defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#ଅନ୍ୟଥା
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ପୂର୍ଣ୍ଣସଂଖ୍ୟା) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#ନଚେତଯଦି defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#ଅନ୍ୟଥା
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  ଅକ୍ଷର trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
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

  ଯଦି (trans == 'N') i = 0;
  ଯଦି (trans == 'T') i = 1;
  ଯଦି (trans == 'R') i = 0;
  ଯଦି (trans == 'C') i = 1;

  ଯଦି (incy == 0)	 info = 13;
  ଯଦି (incx == 0)	 info = 10;
  ଯଦି (lda < kl + ku + 1) info = 8;
  ଯଦି (kl < 0)		 info = 5;
  ଯଦି (ku < 0)		 info = 4;
  ଯଦି (n < 0)		 info = 3;
  ଯଦି (m < 0)		 info = 2;
  ଯଦି (i < 0)		 info = 1;

  trans = i;

  ଯଦି (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order,
	   ଗଣନା CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (TransA == CblasNoTrans)     trans = 0;
    ଯଦି (TransA == CblasTrans)       trans = 1;
    ଯଦି (TransA == CblasConjNoTrans) trans = 0;
    ଯଦି (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ଯଦି (incy == 0)	 info = 13;
    ଯଦି (incx == 0)	 info = 10;
    ଯଦି (lda < kl + ku + 1) info = 8;
    ଯଦି (kl < 0)		 info = 5;
    ଯଦି (ku < 0)		 info = 4;
    ଯଦି (n < 0)		 info = 3;
    ଯଦି (m < 0)		 info = 2;
    ଯଦି (trans < 0)	 info = 1;
  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (TransA == CblasNoTrans)     trans = 1;
    ଯଦି (TransA == CblasTrans)       trans = 0;
    ଯଦି (TransA == CblasConjNoTrans) trans = 1;
    ଯଦି (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ଯଦି (incy == 0)	 info = 13;
    ଯଦି (incx == 0)	 info = 10;
    ଯଦି (lda < kl + ku + 1) info = 8;
    ଯଦି (kl < 0)		 info = 5;
    ଯଦି (ku < 0)		 info = 4;
    ଯଦି (n < 0)		 info = 3;
    ଯଦି (m < 0)		 info = 2;
    ଯଦି (trans < 0)	 info = 1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#endif

  ଯଦି ((m==0) || (n==0)) ଫେରନ୍ତୁ;

  lenx = n;
  leny = m;
  ଯଦି (trans) lenx = m;
  ଯଦି (trans) leny = n;

  ଯଦି (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ଖାଲି, 0, ଖାଲି, 0);

  ଯଦି (alpha == ZERO) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0) x -= (lenx-1)*incx;
  ଯଦି (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ଯଦି (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ifdef SMP
  ଯଦି (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  ଅନ୍ୟଥା
    nthreads = ଅଙ୍କ_cpu_avail(2);

  ଯଦି (nthreads == 1) {
#endif

  (gbmv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ଅନ୍ୟଥା {

    (gbmv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
