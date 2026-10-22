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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#अञ्ञथा
#define ERROR_NAME "SGBMV "
#endif

स्थिर शून्य (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#नोचेद्यदि defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#अञ्ञथा
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#नोचेद्यदि defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#अञ्ञथा
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

शून्य NAME(वर्ण *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  वर्ण trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
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

  सचे (trans == 'N') i = 0;
  सचे (trans == 'T') i = 1;
  सचे (trans == 'R') i = 0;
  सचे (trans == 'C') i = 1;

  सचे (incy == 0)	 info = 13;
  सचे (incx == 0)	 info = 10;
  सचे (lda < kl + ku + 1) info = 8;
  सचे (kl < 0)		 info = 5;
  सचे (ku < 0)		 info = 4;
  सचे (n < 0)		 info = 3;
  सचे (m < 0)		 info = 2;
  सचे (i < 0)		 info = 1;

  trans = i;

  सचे (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  पूर्णाङ्क trans;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  सचे (order == CblasColMajor) {
    सचे (TransA == CblasNoTrans)     trans = 0;
    सचे (TransA == CblasTrans)       trans = 1;
    सचे (TransA == CblasConjNoTrans) trans = 0;
    सचे (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    सचे (incy == 0)	 info = 13;
    सचे (incx == 0)	 info = 10;
    सचे (lda < kl + ku + 1) info = 8;
    सचे (kl < 0)		 info = 5;
    सचे (ku < 0)		 info = 4;
    सचे (n < 0)		 info = 3;
    सचे (m < 0)		 info = 2;
    सचे (trans < 0)	 info = 1;
  }

  सचे (order == CblasRowMajor) {
    सचे (TransA == CblasNoTrans)     trans = 1;
    सचे (TransA == CblasTrans)       trans = 0;
    सचे (TransA == CblasConjNoTrans) trans = 1;
    सचे (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    सचे (incy == 0)	 info = 13;
    सचे (incx == 0)	 info = 10;
    सचे (lda < kl + ku + 1) info = 8;
    सचे (kl < 0)		 info = 5;
    सचे (ku < 0)		 info = 4;
    सचे (n < 0)		 info = 3;
    सचे (m < 0)		 info = 2;
    सचे (trans < 0)	 info = 1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  सचे ((m==0) || (n==0)) पटिनिवत्त;

  lenx = n;
  leny = m;
  सचे (trans) lenx = m;
  सचे (trans) leny = n;

  सचे (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  सचे (alpha == ZERO) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0) x -= (lenx-1)*incx;
  सचे (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  सचे (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#ifdef SMP
  सचे (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  अञ्ञथा
    nthreads = अङ्कः_cpu_avail(2);

  सचे (nthreads == 1) {
#endif

  (gbmv[(पूर्णाङ्क)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } अञ्ञथा {

    (gbmv_thread[(पूर्णाङ्क)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  पटिनिवत्त;
}
