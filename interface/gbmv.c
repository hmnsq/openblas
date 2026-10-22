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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#இல்லையெனில்
#define ERROR_NAME "SGBMV "
#endif

நிலையான வெறுமை (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, வெறுமை *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#அல்லதெனில் defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#இல்லையெனில்
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
நிலையான முழு (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, முழு) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#அல்லதெனில் defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#இல்லையெனில்
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

வெறுமை NAME(எழுத்து *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  எழுத்து trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  முழு nthreads;
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

  எனில் (trans == 'N') i = 0;
  எனில் (trans == 'T') i = 1;
  எனில் (trans == 'R') i = 0;
  எனில் (trans == 'C') i = 1;

  எனில் (incy == 0)	 info = 13;
  எனில் (incx == 0)	 info = 10;
  எனில் (lda < kl + ku + 1) info = 8;
  எனில் (kl < 0)		 info = 5;
  எனில் (ku < 0)		 info = 4;
  எனில் (n < 0)		 info = 3;
  எனில் (m < 0)		 info = 2;
  எனில் (i < 0)		 info = 1;

  trans = i;

  எனில் (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order,
	   எண்ணிக்கை CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  முழு trans;
#ifdef SMP
  முழு nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (TransA == CblasNoTrans)     trans = 0;
    எனில் (TransA == CblasTrans)       trans = 1;
    எனில் (TransA == CblasConjNoTrans) trans = 0;
    எனில் (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    எனில் (incy == 0)	 info = 13;
    எனில் (incx == 0)	 info = 10;
    எனில் (lda < kl + ku + 1) info = 8;
    எனில் (kl < 0)		 info = 5;
    எனில் (ku < 0)		 info = 4;
    எனில் (n < 0)		 info = 3;
    எனில் (m < 0)		 info = 2;
    எனில் (trans < 0)	 info = 1;
  }

  எனில் (order == CblasRowMajor) {
    எனில் (TransA == CblasNoTrans)     trans = 1;
    எனில் (TransA == CblasTrans)       trans = 0;
    எனில் (TransA == CblasConjNoTrans) trans = 1;
    எனில் (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    எனில் (incy == 0)	 info = 13;
    எனில் (incx == 0)	 info = 10;
    எனில் (lda < kl + ku + 1) info = 8;
    எனில் (kl < 0)		 info = 5;
    எனில் (ku < 0)		 info = 4;
    எனில் (n < 0)		 info = 3;
    எனில் (m < 0)		 info = 2;
    எனில் (trans < 0)	 info = 1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#endif

  எனில் ((m==0) || (n==0)) திரும்பு;

  lenx = n;
  leny = m;
  எனில் (trans) lenx = m;
  எனில் (trans) leny = n;

  எனில் (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), வெற்று, 0, வெற்று, 0);

  எனில் (alpha == ZERO) திரும்பு;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  எனில் (incx < 0) x -= (lenx-1)*incx;
  எனில் (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  எனில் (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#ifdef SMP
  எனில் (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  இல்லையெனில்
    nthreads = இலக்கம்_cpu_avail(2);

  எனில் (nthreads == 1) {
#endif

  (gbmv[(முழு)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } இல்லையெனில் {

    (gbmv_thread[(முழு)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  திரும்பு;
}
