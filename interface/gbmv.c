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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#નહીંતર
#define ERROR_NAME "SGBMV "
#endif

સ્થિર શૂન્ય (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, શૂન્ય *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#અથવાજો defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#નહીંતર
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
સ્થિર પૂર્ણાંક (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, પૂર્ણાંક) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#અથવાજો defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#નહીંતર
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

શૂન્ય NAME(અક્ષર *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  અક્ષર trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  પૂર્ણાંક nthreads;
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

  જો (trans == 'N') i = 0;
  જો (trans == 'T') i = 1;
  જો (trans == 'R') i = 0;
  જો (trans == 'C') i = 1;

  જો (incy == 0)	 info = 13;
  જો (incx == 0)	 info = 10;
  જો (lda < kl + ku + 1) info = 8;
  જો (kl < 0)		 info = 5;
  જો (ku < 0)		 info = 4;
  જો (n < 0)		 info = 3;
  જો (m < 0)		 info = 2;
  જો (i < 0)		 info = 1;

  trans = i;

  જો (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order,
	   ગણના CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  પૂર્ણાંક trans;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  જો (order == CblasColMajor) {
    જો (TransA == CblasNoTrans)     trans = 0;
    જો (TransA == CblasTrans)       trans = 1;
    જો (TransA == CblasConjNoTrans) trans = 0;
    જો (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    જો (incy == 0)	 info = 13;
    જો (incx == 0)	 info = 10;
    જો (lda < kl + ku + 1) info = 8;
    જો (kl < 0)		 info = 5;
    જો (ku < 0)		 info = 4;
    જો (n < 0)		 info = 3;
    જો (m < 0)		 info = 2;
    જો (trans < 0)	 info = 1;
  }

  જો (order == CblasRowMajor) {
    જો (TransA == CblasNoTrans)     trans = 1;
    જો (TransA == CblasTrans)       trans = 0;
    જો (TransA == CblasConjNoTrans) trans = 1;
    જો (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    જો (incy == 0)	 info = 13;
    જો (incx == 0)	 info = 10;
    જો (lda < kl + ku + 1) info = 8;
    જો (kl < 0)		 info = 5;
    જો (ku < 0)		 info = 4;
    જો (n < 0)		 info = 3;
    જો (m < 0)		 info = 2;
    જો (trans < 0)	 info = 1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  જો ((m==0) || (n==0)) પાછા;

  lenx = n;
  leny = m;
  જો (trans) lenx = m;
  જો (trans) leny = n;

  જો (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ખાલી, 0, ખાલી, 0);

  જો (alpha == ZERO) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0) x -= (lenx-1)*incx;
  જો (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#ifdef SMP
  જો (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  નહીંતર
    nthreads = અંક_cpu_avail(2);

  જો (nthreads == 1) {
#endif

  (gbmv[(પૂર્ણાંક)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } નહીંતર {

    (gbmv_thread[(પૂર્ણાંક)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  પાછા;
}
