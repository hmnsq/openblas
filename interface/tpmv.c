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
#define ERROR_NAME "QTPMV "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DTPMV "
#അല്ലെങ്കിൽ
#define ERROR_NAME "STPMV "
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, ശൂന്യം *) = {
#ifdef XDOUBLE
  qtpmv_NUU, qtpmv_NUN, qtpmv_NLU, qtpmv_NLN,
  qtpmv_TUU, qtpmv_TUN, qtpmv_TLU, qtpmv_TLN,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  dtpmv_NUU, dtpmv_NUN, dtpmv_NLU, dtpmv_NLN,
  dtpmv_TUU, dtpmv_TUN, dtpmv_TLU, dtpmv_TLN,
#അല്ലെങ്കിൽ
  stpmv_NUU, stpmv_NUN, stpmv_NLU, stpmv_NLN,
  stpmv_TUU, stpmv_TUN, stpmv_TLU, stpmv_TLN,
#endif
};

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
  qtpmv_thread_NUU, qtpmv_thread_NUN, qtpmv_thread_NLU, qtpmv_thread_NLN,
  qtpmv_thread_TUU, qtpmv_thread_TUN, qtpmv_thread_TLU, qtpmv_thread_TLN,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  dtpmv_thread_NUU, dtpmv_thread_NUN, dtpmv_thread_NLU, dtpmv_thread_NLN,
  dtpmv_thread_TUU, dtpmv_thread_TUN, dtpmv_thread_TLU, dtpmv_thread_TLN,
#അല്ലെങ്കിൽ
  stpmv_thread_NUU, stpmv_thread_NUN, stpmv_thread_NLU, stpmv_thread_NLN,
  stpmv_thread_TUU, stpmv_thread_TUN, stpmv_thread_TLU, stpmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, അക്ഷരം *TRANS, അക്ഷരം *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  അക്ഷരം uplo_arg  = *UPLO;
  അക്ഷരം trans_arg = *TRANS;
  അക്ഷരം diag_arg  = *DIAG;

  blasint n    = *N;
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

  എങ്കിൽ (incx == 0)          info =  7;
  എങ്കിൽ (n < 0)              info =  4;
  എങ്കിൽ (unit  < 0)          info =  3;
  എങ്കിൽ (trans < 0)          info =  2;
  എങ്കിൽ (uplo  < 0)          info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ
#ifndef COMPLEX
ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA, എണ്ണൽ CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {
#അല്ലെങ്കിൽ
ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA, എണ്ണൽ CBLAS_DIAG Diag,
	   blasint n, ശൂന്യം  *va, ശൂന്യം  *vx, blasint incx) {
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
#endif

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

    എങ്കിൽ (incx == 0)          info =  7;
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

    എങ്കിൽ (incx == 0)          info =  7;
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

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

#ifdef SMP
  } അല്ലെങ്കിൽ {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  തിരികെ;
}
