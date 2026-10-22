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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#ꦭꦶꦪ
#define ERROR_NAME "SGBMV "
#endif

static void (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#elif defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#ꦭꦶꦪ
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
static puur_naa_nka (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#elif defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#ꦭꦶꦪ
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

void NAME(var_na *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  var_na trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  puur_naa_nka nthreads;
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

  ꦪꦼꦤ꧀ (trans == 'N') i = 0;
  ꦪꦼꦤ꧀ (trans == 'T') i = 1;
  ꦪꦼꦤ꧀ (trans == 'R') i = 0;
  ꦪꦼꦤ꧀ (trans == 'C') i = 1;

  ꦪꦼꦤ꧀ (incy == 0)	 info = 13;
  ꦪꦼꦤ꧀ (incx == 0)	 info = 10;
  ꦪꦼꦤ꧀ (lda < kl + ku + 1) info = 8;
  ꦪꦼꦤ꧀ (kl < 0)		 info = 5;
  ꦪꦼꦤ꧀ (ku < 0)		 info = 4;
  ꦪꦼꦤ꧀ (n < 0)		 info = 3;
  ꦪꦼꦤ꧀ (m < 0)		 info = 2;
  ꦪꦼꦤ꧀ (i < 0)		 info = 1;

  trans = i;

  ꦪꦼꦤ꧀ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ꦮꦭꦶ;
  }

#ꦭꦶꦪ

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  puur_naa_nka trans;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ꦪꦼꦤ꧀ (order == CblasColMajor) {
    ꦪꦼꦤ꧀ (TransA == CblasNoTrans)     trans = 0;
    ꦪꦼꦤ꧀ (TransA == CblasTrans)       trans = 1;
    ꦪꦼꦤ꧀ (TransA == CblasConjNoTrans) trans = 0;
    ꦪꦼꦤ꧀ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ꦪꦼꦤ꧀ (incy == 0)	 info = 13;
    ꦪꦼꦤ꧀ (incx == 0)	 info = 10;
    ꦪꦼꦤ꧀ (lda < kl + ku + 1) info = 8;
    ꦪꦼꦤ꧀ (kl < 0)		 info = 5;
    ꦪꦼꦤ꧀ (ku < 0)		 info = 4;
    ꦪꦼꦤ꧀ (n < 0)		 info = 3;
    ꦪꦼꦤ꧀ (m < 0)		 info = 2;
    ꦪꦼꦤ꧀ (trans < 0)	 info = 1;
  }

  ꦪꦼꦤ꧀ (order == CblasRowMajor) {
    ꦪꦼꦤ꧀ (TransA == CblasNoTrans)     trans = 1;
    ꦪꦼꦤ꧀ (TransA == CblasTrans)       trans = 0;
    ꦪꦼꦤ꧀ (TransA == CblasConjNoTrans) trans = 1;
    ꦪꦼꦤ꧀ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    ꦪꦼꦤ꧀ (incy == 0)	 info = 13;
    ꦪꦼꦤ꧀ (incx == 0)	 info = 10;
    ꦪꦼꦤ꧀ (lda < kl + ku + 1) info = 8;
    ꦪꦼꦤ꧀ (kl < 0)		 info = 5;
    ꦪꦼꦤ꧀ (ku < 0)		 info = 4;
    ꦪꦼꦤ꧀ (n < 0)		 info = 3;
    ꦪꦼꦤ꧀ (m < 0)		 info = 2;
    ꦪꦼꦤ꧀ (trans < 0)	 info = 1;
  }

  ꦪꦼꦤ꧀ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ꦮꦭꦶ;
  }

#endif

  ꦪꦼꦤ꧀ ((m==0) || (n==0)) ꦮꦭꦶ;

  lenx = n;
  leny = m;
  ꦪꦼꦤ꧀ (trans) lenx = m;
  ꦪꦼꦤ꧀ (trans) leny = n;

  ꦪꦼꦤ꧀ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  ꦪꦼꦤ꧀ (alpha == ZERO) ꦮꦭꦶ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ꦪꦼꦤ꧀ (incx < 0) x -= (lenx-1)*incx;
  ꦪꦼꦤ꧀ (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ꦪꦼꦤ꧀ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ꦮꦭꦶ;
  }

#ifdef SMP
  ꦪꦼꦤ꧀ (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  ꦭꦶꦪ
    nthreads = num_cpu_avail(2);

  ꦪꦼꦤ꧀ (nthreads == 1) {
#endif

  (gbmv[(puur_naa_nka)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ꦭꦶꦪ {

    (gbmv_thread[(puur_naa_nka)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  ꦮꦭꦶ;
}
