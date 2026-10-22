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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGBMV "
#endif

സ്ഥിതം ശൂന്യം (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ശൂന്യം *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#അല്ലെങ്കിൽ
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#അല്ലെങ്കിൽ
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  അക്ഷരം trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
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

  എങ്കിൽ (trans == 'N') i = 0;
  എങ്കിൽ (trans == 'T') i = 1;
  എങ്കിൽ (trans == 'R') i = 0;
  എങ്കിൽ (trans == 'C') i = 1;

  എങ്കിൽ (incy == 0)	 info = 13;
  എങ്കിൽ (incx == 0)	 info = 10;
  എങ്കിൽ (lda < kl + ku + 1) info = 8;
  എങ്കിൽ (kl < 0)		 info = 5;
  എങ്കിൽ (ku < 0)		 info = 4;
  എങ്കിൽ (n < 0)		 info = 3;
  എങ്കിൽ (m < 0)		 info = 2;
  എങ്കിൽ (i < 0)		 info = 1;

  trans = i;

  എങ്കിൽ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  പൂർണ്ണസംഖ്യ trans;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (TransA == CblasNoTrans)     trans = 0;
    എങ്കിൽ (TransA == CblasTrans)       trans = 1;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 0;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    എങ്കിൽ (incy == 0)	 info = 13;
    എങ്കിൽ (incx == 0)	 info = 10;
    എങ്കിൽ (lda < kl + ku + 1) info = 8;
    എങ്കിൽ (kl < 0)		 info = 5;
    എങ്കിൽ (ku < 0)		 info = 4;
    എങ്കിൽ (n < 0)		 info = 3;
    എങ്കിൽ (m < 0)		 info = 2;
    എങ്കിൽ (trans < 0)	 info = 1;
  }

  എങ്കിൽ (order == CblasRowMajor) {
    എങ്കിൽ (TransA == CblasNoTrans)     trans = 1;
    എങ്കിൽ (TransA == CblasTrans)       trans = 0;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 1;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    എങ്കിൽ (incy == 0)	 info = 13;
    എങ്കിൽ (incx == 0)	 info = 10;
    എങ്കിൽ (lda < kl + ku + 1) info = 8;
    എങ്കിൽ (kl < 0)		 info = 5;
    എങ്കിൽ (ku < 0)		 info = 4;
    എങ്കിൽ (n < 0)		 info = 3;
    എങ്കിൽ (m < 0)		 info = 2;
    എങ്കിൽ (trans < 0)	 info = 1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#endif

  എങ്കിൽ ((m==0) || (n==0)) തിരികെ;

  lenx = n;
  leny = m;
  എങ്കിൽ (trans) lenx = m;
  എങ്കിൽ (trans) leny = n;

  എങ്കിൽ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ശൂന്യ, 0, ശൂന്യ, 0);

  എങ്കിൽ (alpha == ZERO) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0) x -= (lenx-1)*incx;
  എങ്കിൽ (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  എങ്കിൽ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#ifdef SMP
  എങ്കിൽ (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  അല്ലെങ്കിൽ
    nthreads = അക്കം_cpu_avail(2);

  എങ്കിൽ (nthreads == 1) {
#endif

  (gbmv[(പൂർണ്ണസംഖ്യ)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } അല്ലെങ്കിൽ {

    (gbmv_thread[(പൂർണ്ണസംഖ്യ)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  തിരികെ;
}
