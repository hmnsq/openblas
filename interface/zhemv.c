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
#иначе
#define ERROR_NAME "CHEMV "
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  char uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;
#ifdef SMP
  целое nthreads;
#endif

  целое (*hemv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    HEMV_U, HEMV_L, HEMV_V, HEMV_M,
  };

#ifdef SMP
  целое (*hemv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, целое) = {
    HEMV_THREAD_U, HEMV_THREAD_L, HEMV_THREAD_V, HEMV_THREAD_M,
  };
#endif

  blasint info;
  целое uplo;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  если (uplo_arg  == 'U') uplo  = 0;
  если (uplo_arg  == 'L') uplo  = 1;
  если (uplo_arg  == 'V') uplo  = 2;
  если (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  если (incy == 0)          info = 10;
  если (incx == 0)          info =  7;
  если (lda  < MAX(1, n))   info =  5;
  если (n < 0)              info =  2;
  если (uplo  < 0)          info =  1;

  если (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#иначе

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint n, void *VALPHA,
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
  целое uplo;
  blasint info;
#ifdef SMP
  целое nthreads;
#endif

  целое (*hemv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    HEMV_U, HEMV_L, HEMV_V, HEMV_M,
  };

#ifdef SMP
  целое (*hemv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, целое) = {
    HEMV_THREAD_U, HEMV_THREAD_L, HEMV_THREAD_V, HEMV_THREAD_M,
  };
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  если (order == CblasColMajor) {

    если (Uplo == CblasUpper) uplo  = 0;
    если (Uplo == CblasLower) uplo  = 1;

    info = -1;

    если (incy == 0)          info = 10;
    если (incx == 0)          info =  7;
    если (lda  < MAX(1, n))   info =  5;
    если (n < 0)              info =  2;
    если (uplo  < 0)          info =  1;
  }

  если (order == CblasRowMajor) {

    если (Uplo == CblasUpper) uplo  = 3;
    если (Uplo == CblasLower) uplo  = 2;

    info = -1;

    если (incy == 0)          info = 10;
    если (incx == 0)          info =  7;
    если (lda  < MAX(1, n))   info =  5;
    если (n < 0)              info =  2;
    если (uplo  < 0)          info =  1;
  }

  если (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#endif

  если (n == 0) вернуть;

  если ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), NULL, 0, NULL, 0);

  если ((alpha_r == ZERO) && (alpha_i == ZERO)) вернуть;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  если (incx < 0 ) x -= (n - 1) * incx * 2;
  если (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  если (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#ifdef SMP
  если (n<MULTI_THREAD_MINIMAL) {
	  nthreads = 1 ;
  } иначе {
  	  nthreads = num_cpu_avail(2);
  };

  если (nthreads == 1) {
#endif

  (hemv[uplo])(n, n, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } иначе {

    (hemv_thread[uplo])(n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, 2 * n * n);

  IDEBUG_END;

  вернуть;
}
