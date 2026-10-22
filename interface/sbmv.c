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
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ifdef XDOUBLE
#define ERROR_NAME "QSBMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DSBMV "
#ᮜᮤᮚ
#define ERROR_NAME "SSBMV "
#endif

static  puur_naa_nka (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qsbmv_U, qsbmv_L,
#elif defined(DOUBLE)
  dsbmv_U, dsbmv_L,
#ᮜᮤᮚ
  ssbmv_U, ssbmv_L,
#endif
};

#ifdef SMPTEST
static  puur_naa_nka (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  qsbmv_thread_U, qsbmv_thread_L,
#elif defined(DOUBLE)
  dsbmv_thread_U, dsbmv_thread_L,
#ᮜᮤᮚ
  ssbmv_thread_U, ssbmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  var_na uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  puur_naa_nka uplo;
  FLOAT *buffer;
#ifdef SMPTEST
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ᮚᮦᮔ᮪ (uplo_arg  == 'U') uplo  = 0;
  ᮚᮦᮔ᮪ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᮚᮦᮔ᮪ (incy == 0)          info = 11;
  ᮚᮦᮔ᮪ (incx == 0)          info =  8;
  ᮚᮦᮔ᮪ (lda  < k + 1)       info =  6;
  ᮚᮦᮔ᮪ (k < 0)              info =  3;
  ᮚᮦᮔ᮪ (n < 0)              info =  2;
  ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;

  ᮚᮦᮔ᮪ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#ᮜᮤᮚ

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  puur_naa_nka uplo;
  blasint info;
#ifdef SMPTEST
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ᮚᮦᮔ᮪ (order == CblasColMajor) {
    ᮚᮦᮔ᮪ (Uplo == CblasUpper)         uplo  = 0;
    ᮚᮦᮔ᮪ (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ᮚᮦᮔ᮪ (incy == 0)          info = 11;
    ᮚᮦᮔ᮪ (incx == 0)          info =  8;
    ᮚᮦᮔ᮪ (lda  < k + 1)       info =  6;
    ᮚᮦᮔ᮪ (k < 0)              info =  3;
    ᮚᮦᮔ᮪ (n < 0)              info =  2;
    ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;
  }

  ᮚᮦᮔ᮪ (order == CblasRowMajor) {
    ᮚᮦᮔ᮪ (Uplo == CblasUpper)         uplo  = 1;
    ᮚᮦᮔ᮪ (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    ᮚᮦᮔ᮪ (incy == 0)          info = 11;
    ᮚᮦᮔ᮪ (incx == 0)          info =  8;
    ᮚᮦᮔ᮪ (lda  < k + 1)       info =  6;
    ᮚᮦᮔ᮪ (k < 0)              info =  3;
    ᮚᮦᮔ᮪ (n < 0)              info =  2;
    ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;
  }

  ᮚᮦᮔ᮪ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#endif

  ᮚᮦᮔ᮪ (n == 0) ᮝᮜᮤ;

  ᮚᮦᮔ᮪ (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  ᮚᮦᮔ᮪ (alpha == ZERO) ᮝᮜᮤ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᮚᮦᮔ᮪ (incx < 0 ) x -= (n - 1) * incx;
  ᮚᮦᮔ᮪ (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᮚᮦᮔ᮪ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#ifdef SMPTEST
  nthreads = num_cpu_avail(2);

  ᮚᮦᮔ᮪ (nthreads == 1) {
#endif

  (sbmv[uplo])(n, k, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMPTEST
  } ᮜᮤᮚ {

    (sbmv_thread[uplo])(n, k, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  ᮝᮜᮤ;
}
