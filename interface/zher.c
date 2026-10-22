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
#define ERROR_NAME "XHER  "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZHER  "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CHER  "
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*her[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xher_U, xher_L, xher_V, xher_M,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  zher_U, zher_L, zher_V, zher_M,
#അല്ലെങ്കിൽ
  cher_U, cher_L, cher_V, cher_M,
#endif
};

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*her_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
  xher_thread_U, xher_thread_L, xher_thread_V, xher_thread_M,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  zher_thread_U, zher_thread_L, zher_thread_V, zher_thread_M,
#അല്ലെങ്കിൽ
  cher_thread_U, cher_thread_L, cher_thread_V, cher_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, blasint *N, FLOAT  *ALPHA,
	  FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  അക്ഷരം uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  പൂർണ്ണസംഖ്യ uplo;
  FLOAT *buffer;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  എങ്കിൽ (uplo_arg  == 'U') uplo  = 0;
  എങ്കിൽ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  എങ്കിൽ (lda  < MAX(1, n))   info =  7;
  എങ്കിൽ (incx == 0)          info =  5;
  എങ്കിൽ (n < 0)              info =  2;
  എങ്കിൽ (uplo  < 0)          info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo, blasint n, FLOAT alpha, ശൂന്യം *vx, blasint incx, ശൂന്യം *va, blasint lda) {

  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;

  FLOAT *buffer;
  പൂർണ്ണസംഖ്യ uplo;
  blasint info;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {

    എങ്കിൽ (Uplo == CblasUpper) uplo  = 0;
    എങ്കിൽ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    എങ്കിൽ (lda  < MAX(1, n))   info =  7;
    എങ്കിൽ (incx == 0)          info =  5;
    എങ്കിൽ (n < 0)              info =  2;
    എങ്കിൽ (uplo  < 0)          info =  1;

  }

  എങ്കിൽ (order == CblasRowMajor) {

    എങ്കിൽ (Uplo == CblasUpper) uplo  = 3;
    എങ്കിൽ (Uplo == CblasLower) uplo  = 2;

    info = -1;

    എങ്കിൽ (lda  < MAX(1, n))   info =  7;
    എങ്കിൽ (incx == 0)          info =  5;
    എങ്കിൽ (n < 0)              info =  2;
    എങ്കിൽ (uplo  < 0)          info =  1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#endif

  എങ്കിൽ (n == 0) തിരികെ;

  എങ്കിൽ (alpha == ZERO) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  എങ്കിൽ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#ifdef SMP
  എങ്കിൽ (n < 100)
    nthreads = 1;
  അല്ലെങ്കിൽ
    nthreads = അക്കം_cpu_avail(2);

  എങ്കിൽ (nthreads == 1) {
#endif

  (her[uplo])(n, alpha, x, incx, a, lda, buffer);

#ifdef SMP
  } അല്ലെങ്കിൽ {

    (her_thread[uplo])(n, alpha, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  തിരികെ;
}
