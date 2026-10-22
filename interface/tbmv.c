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
#അന്തര്ഭാവ <ctype.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTBMV "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DTBMV "
#അല്ലെങ്കിൽ
#define ERROR_NAME "STBMV "
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ശൂന്യം *) = {
#ifdef XDOUBLE
  qtbmv_NUU, qtbmv_NUN, qtbmv_NLU, qtbmv_NLN,
  qtbmv_TUU, qtbmv_TUN, qtbmv_TLU, qtbmv_TLN,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  dtbmv_NUU, dtbmv_NUN, dtbmv_NLU, dtbmv_NLN,
  dtbmv_TUU, dtbmv_TUN, dtbmv_TLU, dtbmv_TLN,
#അല്ലെങ്കിൽ
  stbmv_NUU, stbmv_NUN, stbmv_NLU, stbmv_NLN,
  stbmv_TUU, stbmv_TUN, stbmv_TLU, stbmv_TLN,
#endif
};

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
  qtbmv_thread_NUU, qtbmv_thread_NUN, qtbmv_thread_NLU, qtbmv_thread_NLN,
  qtbmv_thread_TUU, qtbmv_thread_TUN, qtbmv_thread_TLU, qtbmv_thread_TLN,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  dtbmv_thread_NUU, dtbmv_thread_NUN, dtbmv_thread_NLU, dtbmv_thread_NLN,
  dtbmv_thread_TUU, dtbmv_thread_TUN, dtbmv_thread_TLU, dtbmv_thread_TLN,
#അല്ലെങ്കിൽ
  stbmv_thread_NUU, stbmv_thread_NUN, stbmv_thread_NLU, stbmv_thread_NLN,
  stbmv_thread_TUU, stbmv_thread_TUN, stbmv_thread_TLU, stbmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, അക്ഷരം *TRANS, അക്ഷരം *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  അക്ഷരം uplo_arg  = *UPLO;
  അക്ഷരം trans_arg = *TRANS;
  അക്ഷരം diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  പൂർണ്ണസംഖ്യ uplo;
  പൂർണ്ണസംഖ്യ unit;
  പൂർണ്ണസംഖ്യ trans;
  FLOAT *buffer;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  എങ്കിൽ (trans_arg == 'N') trans = 0;
  എങ്കിൽ (trans_arg == 'T') trans = 1;
  എങ്കിൽ (trans_arg == 'R') trans = 0;
  എങ്കിൽ (trans_arg == 'C') trans = 1;

  എങ്കിൽ (diag_arg  == 'U') unit  = 0;
  എങ്കിൽ (diag_arg  == 'N') unit  = 1;

  എങ്കിൽ (uplo_arg  == 'U') uplo  = 0;
  എങ്കിൽ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  എങ്കിൽ (incx == 0)          info =  9;
  എങ്കിൽ (lda < k + 1)        info =  7;
  എങ്കിൽ (k < 0)              info =  5;
  എങ്കിൽ (n < 0)              info =  4;
  എങ്കിൽ (unit  < 0)          info =  3;
  എങ്കിൽ (trans < 0)          info =  2;
  എങ്കിൽ (uplo  < 0)          info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA, എണ്ണൽ CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  പൂർണ്ണസംഖ്യ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (Uplo == CblasUpper)         uplo  = 0;
    എങ്കിൽ (Uplo == CblasLower)         uplo  = 1;

    എങ്കിൽ (TransA == CblasNoTrans)     trans = 0;
    എങ്കിൽ (TransA == CblasTrans)       trans = 1;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 0;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 1;

    എങ്കിൽ (Diag == CblasUnit)          unit  = 0;
    എങ്കിൽ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    എങ്കിൽ (incx == 0)          info =  9;
    എങ്കിൽ (lda < k + 1)        info =  7;
    എങ്കിൽ (k < 0)              info =  5;
    എങ്കിൽ (n < 0)              info =  4;
    എങ്കിൽ (unit  < 0)          info =  3;
    എങ്കിൽ (trans < 0)          info =  2;
    എങ്കിൽ (uplo  < 0)          info =  1;
  }

  എങ്കിൽ (order == CblasRowMajor) {
    എങ്കിൽ (Uplo == CblasUpper)         uplo  = 1;
    എങ്കിൽ (Uplo == CblasLower)         uplo  = 0;

    എങ്കിൽ (TransA == CblasNoTrans)     trans = 1;
    എങ്കിൽ (TransA == CblasTrans)       trans = 0;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 1;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 0;

    എങ്കിൽ (Diag == CblasUnit)          unit  = 0;
    എങ്കിൽ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    എങ്കിൽ (incx == 0)          info =  9;
    എങ്കിൽ (lda < k + 1)        info =  7;
    എങ്കിൽ (k < 0)              info =  5;
    എങ്കിൽ (n < 0)              info =  4;
    എങ്കിൽ (unit  < 0)          info =  3;
    എങ്കിൽ (trans < 0)          info =  2;
    എങ്കിൽ (uplo  < 0)          info =  1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#endif

  എങ്കിൽ (n == 0) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  എങ്കിൽ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#ifdef SMP
  nthreads = അക്കം_cpu_avail(2);

  എങ്കിൽ (n < 50) nthreads = 1;
  എങ്കിൽ (nthreads > 2 && n < 500) nthreads = 2;

  എങ്കിൽ (nthreads == 1) {
#endif

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

#ifdef SMP
  } അല്ലെങ്കിൽ {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  തിരികെ;
}
