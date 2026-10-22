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

// this is smallest dimension N of square input a to permit threading
// see graph in issue #1820 for explanation
#define MULTI_THREAD_MINIMAL 362

#ifdef XDOUBLE
#define ERROR_NAME "XHEMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHEMV "
#ᮜᮤᮚ
#define ERROR_NAME "CHEMV "
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  puur_naa_nka (*hemv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    HEMV_U, HEMV_L, HEMV_V, HEMV_M,
  };

#ifdef SMP
  puur_naa_nka (*hemv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
    HEMV_THREAD_U, HEMV_THREAD_L, HEMV_THREAD_V, HEMV_THREAD_M,
  };
#endif

  blasint info;
  puur_naa_nka uplo;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ᮚᮦᮔ᮪ (uplo_arg  == 'U') uplo  = 0;
  ᮚᮦᮔ᮪ (uplo_arg  == 'L') uplo  = 1;
  ᮚᮦᮔ᮪ (uplo_arg  == 'V') uplo  = 2;
  ᮚᮦᮔ᮪ (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  ᮚᮦᮔ᮪ (incy == 0)          info = 10;
  ᮚᮦᮔ᮪ (incx == 0)          info =  7;
  ᮚᮦᮔ᮪ (lda  < MAX(1, n))   info =  5;
  ᮚᮦᮔ᮪ (n < 0)              info =  2;
  ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;

  ᮚᮦᮔ᮪ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#ᮜᮤᮚ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, void *VALPHA,
	   void *va, blasint lda, void *vx, blasint incx, void *VBETA, void *vy, blasint incy) {

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];

  FLOAT *buffer;
  puur_naa_nka uplo;
  blasint info;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  puur_naa_nka (*hemv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    HEMV_U, HEMV_L, HEMV_V, HEMV_M,
  };

#ifdef SMP
  puur_naa_nka (*hemv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
    HEMV_THREAD_U, HEMV_THREAD_L, HEMV_THREAD_V, HEMV_THREAD_M,
  };
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ᮚᮦᮔ᮪ (order == CblasColMajor) {

    ᮚᮦᮔ᮪ (Uplo == CblasUpper) uplo  = 0;
    ᮚᮦᮔ᮪ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ᮚᮦᮔ᮪ (incy == 0)          info = 10;
    ᮚᮦᮔ᮪ (incx == 0)          info =  7;
    ᮚᮦᮔ᮪ (lda  < MAX(1, n))   info =  5;
    ᮚᮦᮔ᮪ (n < 0)              info =  2;
    ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;
  }

  ᮚᮦᮔ᮪ (order == CblasRowMajor) {

    ᮚᮦᮔ᮪ (Uplo == CblasUpper) uplo  = 3;
    ᮚᮦᮔ᮪ (Uplo == CblasLower) uplo  = 2;

    info = -1;

    ᮚᮦᮔ᮪ (incy == 0)          info = 10;
    ᮚᮦᮔ᮪ (incx == 0)          info =  7;
    ᮚᮦᮔ᮪ (lda  < MAX(1, n))   info =  5;
    ᮚᮦᮔ᮪ (n < 0)              info =  2;
    ᮚᮦᮔ᮪ (uplo  < 0)          info =  1;
  }

  ᮚᮦᮔ᮪ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#endif

  ᮚᮦᮔ᮪ (n == 0) ᮝᮜᮤ;

  ᮚᮦᮔ᮪ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), NULL, 0, NULL, 0);

  ᮚᮦᮔ᮪ ((alpha_r == ZERO) && (alpha_i == ZERO)) ᮝᮜᮤ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᮚᮦᮔ᮪ (incx < 0 ) x -= (n - 1) * incx * 2;
  ᮚᮦᮔ᮪ (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᮚᮦᮔ᮪ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᮝᮜᮤ;
  }

#ifdef SMP
  ᮚᮦᮔ᮪ (n<MULTI_THREAD_MINIMAL) {
	  nthreads = 1 ;
  } ᮜᮤᮚ {
  	  nthreads = num_cpu_avail(2);
  };

  ᮚᮦᮔ᮪ (nthreads == 1) {
#endif

  (hemv[uplo])(n, n, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ᮜᮤᮚ {

    (hemv_thread[uplo])(n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, 2 * n * n);

  IDEBUG_END;

  ᮝᮜᮤ;
}
